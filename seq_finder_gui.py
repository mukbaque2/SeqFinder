"""
시퀀싱 결과 서열 검색기 (GUI)
- 찾을 서열을 입력(또는 파일에서 불러오기)
- 시퀀싱 결과 파일(.ab1/.txt/.seq/.fasta)을 여러 개 추가
- 분석 시작 -> 파일별 일치 여부, (ab1 이면) 크로마토그램 품질, 정렬 결과 확인, CSV 저장
"""
import csv
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

from abi_reader import GOOD_LEN, GOOD_Q, Trace, suspicious_reason, weakest_in_run
from analyze_seq import analyze_one, clean_seq, read_seq
from report_image import render_report, render_summary

SEQ_FILETYPES = [("시퀀싱 파일", "*.ab1 *.abi *.txt *.seq *.fasta *.fa *.fas"),
                 ("ab1 원본", "*.ab1 *.abi"), ("모든 파일", "*.*")]
QUERY_FILETYPES = [("서열 파일", "*.txt *.seq *.fasta *.fa *.fas"), ("모든 파일", "*.*")]
ROW_COLORS = {"match": "#d4edda", "diff": "#fff3cd", "none": "#f8d7da", "error": "#e2e3e5"}
BASE_COLORS = {"A": "#1a9850", "C": "#2b6cd4", "G": "#222222", "T": "#d7301f", "N": "#888888"}
FLANK = 12   # 크로마토그램에서 찾은 서열 앞뒤로 더 보여줄 염기 수


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("시퀀싱 서열 검색기")
        self.geometry("1100x840")
        self.minsize(820, 680)
        self.files = []      # 추가된 파일 경로
        self.results = {}    # 파일 경로 -> 분석 결과 dict
        self.last_query = ""

        mono = ("Consolas", 10)
        pad = {"padx": 8, "pady": 3}

        # ---------- 1. 찾을 서열 ----------
        box1 = ttk.LabelFrame(self, text=" 1. 찾을 서열 ")
        box1.pack(fill="x", **pad)
        self.seq_text = tk.Text(box1, height=3, font=mono, wrap="char")
        self.seq_text.pack(side="left", fill="x", expand=True, padx=(8, 4), pady=6)
        self.seq_text.bind("<<Modified>>", self._on_seq_modified)
        side1 = ttk.Frame(box1)
        side1.pack(side="left", fill="y", padx=(0, 8), pady=6)
        ttk.Button(side1, text="파일에서 불러오기", command=self.load_query_file).pack(fill="x")
        ttk.Button(side1, text="지우기", command=lambda: self.seq_text.delete("1.0", "end")).pack(fill="x", pady=4)
        self.seq_len_var = tk.StringVar(value="0 bp")
        ttk.Label(side1, textvariable=self.seq_len_var).pack()

        # ---------- 2. 결과 파일 ----------
        box2 = ttk.LabelFrame(self, text=" 2. 시퀀싱 결과 파일 (.ab1 을 넣으면 크로마토그램 품질도 분석) ")
        box2.pack(fill="x", **pad)
        self.file_list = tk.Listbox(box2, height=4, selectmode="extended")
        self.file_list.pack(side="left", fill="x", expand=True, padx=(8, 4), pady=6)
        side2 = ttk.Frame(box2)
        side2.pack(side="left", fill="y", padx=(0, 8), pady=6)
        ttk.Button(side2, text="파일 추가", command=self.add_files).pack(fill="x")
        ttk.Button(side2, text="선택 삭제", command=self.remove_selected).pack(fill="x", pady=4)
        ttk.Button(side2, text="모두 삭제", command=self.clear_files).pack(fill="x")

        # ---------- 3. 실행 ----------
        bar = ttk.Frame(self)
        bar.pack(fill="x", **pad)
        ttk.Label(bar, text="허용 차이 비율(%)").pack(side="left")
        self.ratio_var = tk.IntVar(value=20)
        ttk.Spinbox(bar, from_=0, to=50, width=5, textvariable=self.ratio_var).pack(side="left", padx=(4, 16))
        self.run_btn = ttk.Button(bar, text="▶ 분석 시작", command=self.run_analysis)
        self.run_btn.pack(side="left")
        self.save_btn = ttk.Button(bar, text="CSV로 저장", command=self.save_csv, state="disabled")
        self.save_btn.pack(side="left", padx=(8, 4))
        self.img_btn = ttk.Button(bar, text="선택 결과 이미지 저장", command=self.save_selected_image, state="disabled")
        self.img_btn.pack(side="left", padx=4)
        self.img_all_btn = ttk.Button(bar, text="전체 이미지 저장", command=self.save_all_images, state="disabled")
        self.img_all_btn.pack(side="left", padx=4)
        self.progress = ttk.Progressbar(bar, mode="determinate", length=160)
        self.progress.pack(side="right")
        self.status_var = tk.StringVar(value="찾을 서열을 입력하고 파일을 추가하세요.")
        ttk.Label(bar, textvariable=self.status_var).pack(side="right", padx=8)

        # ---------- 4. 결과 ----------
        box4 = ttk.LabelFrame(self, text=" 3. 결과 (행을 클릭하면 아래에 정렬과 크로마토그램이 표시됩니다) ")
        box4.pack(fill="both", expand=True, **pad)
        cols = ("file", "verdict", "pos", "site_q", "clean", "grade", "detail")
        self.tree = ttk.Treeview(box4, columns=cols, show="headings", height=6)
        for c, title, w, anchor in [
                ("file", "파일", 150, "w"), ("verdict", "판정", 70, "center"),
                ("pos", "위치", 80, "center"), ("site_q", "찾은 부위 품질", 110, "center"),
                ("clean", "깨끗한 구간", 110, "center"), ("grade", "전체 품질", 70, "center"),
                ("detail", "차이 / 품질 설명", 440, "w")]:
            self.tree.heading(c, text=title)
            self.tree.column(c, width=w, anchor=anchor, stretch=(c == "detail"))
        for tag, color in ROW_COLORS.items():
            self.tree.tag_configure(tag, background=color)
        self.tree.pack(fill="both", expand=True, padx=8, pady=(6, 4))
        self.tree.bind("<<TreeviewSelect>>", self.show_detail)

        self.detail = tk.Text(box4, height=7, font=mono, wrap="none", state="disabled", background="#fafafa")
        self.detail.pack(fill="x", padx=8)

        # 전체 품질 막대 (판독 전체에서 깨끗한 구간 / 찾은 위치)
        ttk.Label(box4, text="전체 품질 (회색 막대 = 염기별 품질, 초록 = 깨끗한 구간, 주황 = 찾은 위치)").pack(anchor="w", padx=8, pady=(6, 0))
        self.overview = tk.Canvas(box4, height=46, background="white", highlightthickness=1,
                                  highlightbackground="#ccc")
        self.overview.pack(fill="x", padx=8)

        # 크로마토그램
        ttk.Label(box4, text="크로마토그램 (찾은 위치 확대, 노란 배경 = 차이 있는 염기, 아래 막대 = 품질 Q, 빨간 선 = Q20)").pack(anchor="w", padx=8, pady=(6, 0))
        chrom_frame = ttk.Frame(box4)
        chrom_frame.pack(fill="x", padx=8, pady=(0, 6))
        self.chrom = tk.Canvas(chrom_frame, height=196, background="white", highlightthickness=1,
                               highlightbackground="#ccc")
        xs = ttk.Scrollbar(chrom_frame, orient="horizontal", command=self.chrom.xview)
        self.chrom.configure(xscrollcommand=xs.set)
        self.chrom.pack(fill="x")
        xs.pack(fill="x")
        self.bind("<Configure>", lambda e: self.after_idle(self.show_detail) if e.widget is self else None)

    # ---------- 입력 ----------
    def _on_seq_modified(self, _event=None):
        self.seq_len_var.set(f"{len(self.get_query())} bp")
        self.seq_text.edit_modified(False)

    def get_query(self):
        return clean_seq(self.seq_text.get("1.0", "end"))

    def load_query_file(self):
        path = filedialog.askopenfilename(title="찾을 서열 파일 선택", filetypes=QUERY_FILETYPES)
        if path:
            self.seq_text.delete("1.0", "end")
            self.seq_text.insert("1.0", read_seq(path))

    def add_files(self):
        paths = filedialog.askopenfilenames(title="시퀀싱 결과 파일 선택", filetypes=SEQ_FILETYPES)
        for p in paths:
            p = Path(p)
            if p not in self.files:
                self.files.append(p)
                self.file_list.insert("end", p.name)
        self.status_var.set(f"파일 {len(self.files)}개")

    def remove_selected(self):
        for i in reversed(self.file_list.curselection()):
            self.file_list.delete(i)
            del self.files[i]
        self.status_var.set(f"파일 {len(self.files)}개")

    def clear_files(self):
        self.file_list.delete(0, "end")
        self.files.clear()
        self.status_var.set("파일 0개")

    # ---------- 분석 ----------
    def run_analysis(self):
        query = self.get_query()
        if not query:
            messagebox.showwarning("입력 필요", "찾을 서열을 입력하세요.")
            return
        bad = sorted(set(query) - set("ACGTN"))
        if bad and not messagebox.askyesno(
                "확인", f"A/C/G/T/N 이외의 문자가 있습니다: {', '.join(bad)}\n그래도 진행할까요?"):
            return
        if not self.files:
            messagebox.showwarning("입력 필요", "시퀀싱 결과 파일을 추가하세요.")
            return
        try:
            ratio = max(0, min(100, int(self.ratio_var.get()))) / 100
        except (tk.TclError, ValueError):
            ratio = 0.2

        self.tree.delete(*self.tree.get_children())
        self.results.clear()
        self._clear_detail()
        self.run_btn.configure(state="disabled")
        for b in (self.save_btn, self.img_btn, self.img_all_btn):
            b.configure(state="disabled")
        self.progress.configure(maximum=len(self.files), value=0)
        self.status_var.set("분석 중...")
        threading.Thread(target=self._worker, args=(query, list(self.files), ratio), daemon=True).start()

    @staticmethod
    def analyze_file(query, f, ratio):
        """파일 하나 분석 -> 결과 dict"""
        trace = None
        if f.suffix.lower() in (".ab1", ".abi"):
            trace = Trace(f)
            read = trace.seq
        else:
            read = read_seq(f)
        r = analyze_one(query, read, ratio)
        r["read_len"] = len(read)
        r["trace"] = trace
        if trace is None:
            return r

        r["qual"] = trace.summary()
        if r["verdict"] != "없음":
            st = trace.region_stats(r["start"], r["end"])
            r["site"] = st
            # 차이 위치마다 품질 근거 (TTT 처럼 같은 염기가 이어지면 그중 가장 약한 염기 기준)
            notes = []
            for label, pos in zip(r["diffs"], r["diff_pos"]):
                p = weakest_in_run(trace, pos)
                notes.append(f"{label}(결과 {p}번: {suspicious_reason(trace, p)})")
            r["diff_notes"] = notes
            r["diff_pos"] = [weakest_in_run(trace, p) for p in r["diff_pos"]]
            if st["min_q"] >= GOOD_Q and not st["mixed"]:
                r["site_label"] = f"좋음 (최소 Q{st['min_q']})"
            else:
                r["site_label"] = f"확인 필요 (최소 Q{st['min_q']})"
        return r

    def _worker(self, query, files, ratio):
        for f in files:
            try:
                r = self.analyze_file(query, f, ratio)
            except Exception as e:   # 파일을 읽을 수 없는 경우 등
                r = {"verdict": "오류", "error": f"{type(e).__name__}: {e}"}
            self.after(0, self._add_result, f, r)
        self.after(0, self._done, query)

    def _add_result(self, f, r):
        self.results[f] = r
        if r["verdict"] == "오류":
            values, tag = (f.name, "오류", "", "", "", "", r["error"]), "error"
        else:
            tag = "match" if r["verdict"] == "일치" else "none" if r["verdict"] == "없음" else "diff"
            found = r["verdict"] != "없음"
            pos = f'{r["start"]}-{r["end"]}' + (" (역)" if r["strand"] == "역상보" else "") if found else ""
            q = r.get("qual")
            clean = f'{q["clean"][0]}-{q["clean"][1]}' if q and q["clean"] else ("-" if not q else "없음")
            grade = q["grade"] if q else "-"
            if not found:
                detail = f'가장 비슷한 곳도 {r["diff"]}개 차이'
            elif r.get("diff_notes"):
                detail = "; ".join(r["diff_notes"])
            else:
                detail = ", ".join(r["diffs"])
            if q:
                detail = (detail + " | " if detail else "") + q["note"]
                if found and r["site"]["mixed"]:
                    detail += f' | 찾은 부위 이중 피크 {len(r["site"]["mixed"])}개'
            site = r.get("site_label", "-") if found else ""
            values = (f.name, r["verdict"], pos, site, clean, grade, detail)
        self.tree.insert("", "end", iid=str(f), values=values, tags=(tag,))
        self.progress.step(1)

    def _done(self, query):
        n = len(self.results)
        ok = sum(r["verdict"] == "일치" for r in self.results.values())
        self.status_var.set(f"완료: {n}개 중 일치 {ok}개 (찾을 서열 {len(query)} bp)")
        self.run_btn.configure(state="normal")
        for b in (self.save_btn, self.img_btn, self.img_all_btn):
            b.configure(state="normal")
        self.last_query = query

    # ---------- 상세 표시 ----------
    def _clear_detail(self):
        self._set_detail("")
        self.overview.delete("all")
        self.chrom.delete("all")

    def show_detail(self, _event=None):
        sel = self.tree.selection()
        if not sel:
            return
        f = Path(sel[0])
        r = self.results.get(f)
        self._clear_detail()
        if not r or r["verdict"] == "오류":
            self._set_detail(r.get("error", "") if r else "")
            return
        found = r["verdict"] != "없음"
        text = f'{f.name}  (결과 서열 {r["read_len"]} bp)\n'
        if found:
            text += (f'판정: {r["verdict"]}  |  방향: {r["strand"]}  |  위치: {r["start"]}-{r["end"]}\n'
                     f'찾을 서열 {r["q_aln"]}\n'
                     f'          {r["bar"]}\n'
                     f'결과 서열 {r["r_aln"]}\n')
            if r.get("diff_notes"):
                text += "차이: " + "; ".join(r["diff_notes"]) + "\n"
            elif r["diffs"]:
                text += "차이: " + ", ".join(r["diffs"]) + "\n"
        else:
            text += f'찾을 서열이 없습니다 (가장 비슷한 곳도 {r["diff"]}개 차이)\n'
        q = r.get("qual")
        if q:
            text += (f'전체 품질: {q["grade"]} | 깨끗한 구간(Q20 기준) {q["clean"][0] if q["clean"] else "-"}'
                     f'-{q["clean"][1] if q["clean"] else "-"} ({q["clean_len"]} bp, 평균 Q{q["clean_meanq"]}) | {q["note"]}\n')
            if found and r["site"]["mixed"]:
                text += "찾은 부위 이중 피크: " + ", ".join(
                    f"{p}번 {b}/{b2} {int(x * 100)}%" for p, b, b2, x in r["site"]["mixed"]) + "\n"
        else:
            text += "(텍스트 파일이라 크로마토그램 품질은 분석하지 않았습니다. .ab1 파일을 넣으면 함께 분석합니다.)\n"
        self._set_detail(text)

        tr = r.get("trace")
        if tr is not None:
            self.draw_overview(tr, q, r if found else None)
            if found:
                self.draw_chromatogram(tr, r)

    def _set_detail(self, text):
        self.detail.configure(state="normal")
        self.detail.delete("1.0", "end")
        self.detail.insert("1.0", text)
        self.detail.configure(state="disabled")

    def draw_overview(self, tr, q, r):
        c = self.overview
        c.update_idletasks()
        w, h = max(c.winfo_width(), 200), int(c["height"])
        n = len(tr)
        sx = (w - 4) / n
        if q and q["clean"]:
            c.create_rectangle(2 + (q["clean"][0] - 1) * sx, 0, 2 + q["clean"][1] * sx, h,
                               fill="#e3f4e6", outline="")
        if r:
            c.create_rectangle(2 + (r["start"] - 1) * sx, 0, 2 + r["end"] * sx, h,
                               fill="#ffd8a8", outline="#f08c00")
        for i, qv in enumerate(tr.qual):
            x = 2 + i * sx
            c.create_line(x, h - 2, x, h - 2 - min(qv, 62) / 62 * (h - 14), fill="#9aa0a6")
        y20 = h - 2 - GOOD_Q / 62 * (h - 14)
        c.create_line(0, y20, w, y20, fill="#e03131", dash=(3, 3))
        for k in range(0, n, 100):
            c.create_text(2 + k * sx, 1, text=str(k or 1), anchor="nw", font=("Arial", 7), fill="#555")

    def draw_chromatogram(self, tr, r):
        c = self.chrom
        c.update_idletasks()
        s = max(1, r["start"] - FLANK)
        e = min(len(tr), r["end"] + FLANK)
        x0, x1 = tr.ploc[s - 1] - 8, tr.ploc[e - 1] + 8
        view_w = max(c.winfo_width(), 300)
        px = max(view_w / (x1 - x0), 1.6)          # 화면 폭에 맞추되 너무 좁으면 스크롤
        total_w = (x1 - x0) * px
        H = int(c["height"])
        top, trace_h, base_y, qbar_top, qbar_h = 6, 112, 130, 142, 48
        X = lambda scan: (scan - x0) * px

        # 차이 위치 / 찾은 구간 배경
        c.create_rectangle(X(tr.ploc[r["start"] - 1]) - 6, 0, X(tr.ploc[r["end"] - 1]) + 6, H,
                           fill="#fff4e6", outline="")
        for p in r.get("diff_pos", []):
            if s <= p <= e:
                xm = X(tr.ploc[p - 1])
                c.create_rectangle(xm - 6 * max(1, px / 1.6), 0, xm + 6 * max(1, px / 1.6), H,
                                   fill="#ffe066", outline="")

        # 신호 크기 (해당 구간 최대값 기준)
        peak = max(max(ch[x0:x1 + 1]) for ch in tr.channels.values()) or 1
        for base, ch in tr.channels.items():
            pts = []
            for scan in range(max(0, x0), min(len(ch), x1 + 1)):
                pts += [X(scan), top + trace_h - ch[scan] / peak * trace_h]
            if len(pts) >= 4:
                c.create_line(*pts, fill=BASE_COLORS.get(base, "#000"), width=1.2)

        # 염기 글자, 번호, 품질 막대
        for i in range(s, e + 1):
            xm = X(tr.ploc[i - 1])
            b = tr.seq[i - 1]
            c.create_text(xm, base_y, text=b, fill=BASE_COLORS.get(b, "#000"), font=("Consolas", 10, "bold"))
            if i % 10 == 0:
                c.create_text(xm, top + 2, text=str(i), anchor="n", font=("Arial", 7), fill="#555")
            qv = tr.qual[i - 1]
            col = "#74c69d" if qv >= GOOD_Q else "#ff8787"
            c.create_rectangle(xm - 4, qbar_top + qbar_h - min(qv, 62) / 62 * qbar_h,
                               xm + 4, qbar_top + qbar_h, fill=col, outline="")
        y20 = qbar_top + qbar_h - GOOD_Q / 62 * qbar_h
        c.create_line(0, y20, total_w, y20, fill="#e03131", dash=(3, 3))
        c.configure(scrollregion=(0, 0, total_w, H))
        c.xview_moveto(0)

    # ---------- 저장 ----------
    def save_csv(self):
        path = filedialog.asksaveasfilename(title="결과 저장", defaultextension=".csv",
                                            initialfile="analysis_result.csv",
                                            filetypes=[("CSV 파일", "*.csv")])
        if not path:
            return
        with open(path, "w", newline="", encoding="utf-8-sig") as fp:   # 엑셀에서 한글 안 깨지게
            w = csv.writer(fp)
            w.writerow(["찾을서열", "파일", "판정", "차이개수", "방향", "시작", "끝", "차이상세",
                        "찾은부위품질", "깨끗한구간시작", "깨끗한구간끝", "깨끗한구간길이",
                        "깨끗한구간평균Q", "전체품질", "품질설명"])
            for f, r in self.results.items():
                if r["verdict"] == "오류":
                    w.writerow([self.last_query, f.name, "오류"] + [""] * 4 + [r["error"]] + [""] * 7)
                    continue
                q = r.get("qual") or {}
                clean = q.get("clean") or ("", "")
                w.writerow([self.last_query, f.name, r["verdict"], r["diff"], r["strand"],
                            r["start"], r["end"], "; ".join(r.get("diff_notes") or r["diffs"]),
                            r.get("site_label", ""), clean[0], clean[1], q.get("clean_len", ""),
                            q.get("clean_meanq", ""), q.get("grade", ""), q.get("note", "")])
        self.status_var.set(f"저장 완료: {Path(path).name}")

    def save_selected_image(self):
        sel = self.tree.selection()
        if not sel:
            messagebox.showinfo("선택 필요", "결과 표에서 저장할 행을 먼저 클릭하세요.")
            return
        f = Path(sel[0])
        path = filedialog.asksaveasfilename(title="결과 이미지 저장", defaultextension=".png",
                                            initialfile=f"{f.stem}_result.png",
                                            filetypes=[("PNG 이미지", "*.png"), ("JPEG 이미지", "*.jpg")])
        if not path:
            return
        try:
            img = render_report(f, self.results[f], self.last_query)
            img.save(path, quality=95) if path.lower().endswith((".jpg", ".jpeg")) else img.save(path)
        except Exception as e:
            messagebox.showerror("저장 실패", f"{type(e).__name__}: {e}")
            return
        self.status_var.set(f"이미지 저장 완료: {Path(path).name}")

    def save_all_images(self):
        folder = filedialog.askdirectory(title="이미지를 저장할 폴더 선택")
        if not folder:
            return
        folder = Path(folder)
        targets = [folder / "00_summary.png"] + [folder / f"{f.stem}_result.png" for f in self.results]
        exists = [p.name for p in targets if p.exists()]
        if exists and not messagebox.askyesno(
                "덮어쓰기 확인", f"이미 있는 파일 {len(exists)}개를 덮어씁니다:\n" + "\n".join(exists[:8])
                + ("\n..." if len(exists) > 8 else "") + "\n\n계속할까요?"):
            return
        try:
            render_summary(self.results, self.last_query).save(targets[0])
            for f, out in zip(self.results, targets[1:]):
                render_report(f, self.results[f], self.last_query).save(out)
                self.status_var.set(f"이미지 저장 중... {out.name}")
                self.update_idletasks()
        except Exception as e:
            messagebox.showerror("저장 실패", f"{type(e).__name__}: {e}")
            return
        self.status_var.set(f"이미지 {len(targets)}개 저장 완료: {folder}")


if __name__ == "__main__":
    App().mainloop()
