"""
분석 결과를 PNG 이미지로 그리기 (Pillow 사용).
- render_report : 파일 하나의 결과 (정렬 + 전체 품질 + 크로마토그램)
- render_summary: 여러 파일의 결과 요약 표
"""
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

from abi_reader import GOOD_Q
from alignment_view import build_lines, combined_aa_lines, combined_dna_lines

W = 1800                 # 이미지 폭(px)
M = 40                   # 좌우 여백
BASES_PER_ROW = 80       # 크로마토그램 한 줄에 그릴 최대 염기 수
FLANK = 12               # 찾은 서열 앞뒤로 더 보여줄 염기 수
BASE_COLORS = {"A": (26, 152, 80), "C": (43, 108, 212), "G": (34, 34, 34), "T": (215, 48, 31), "N": (136, 136, 136)}
VERDICT_BG = {"일치": (212, 237, 218), "없음": (248, 215, 218), "오류": (226, 227, 229)}
DIFF_BG = (255, 243, 205)


def _font(names, size):
    """윈도우 글꼴 중 있는 것을 사용 (한글은 맑은 고딕)"""
    for n in names:
        for p in (Path("C:/Windows/Fonts") / n, Path(n)):
            if p.exists():
                return ImageFont.truetype(str(p), size)
    return ImageFont.load_default()


F_TITLE = _font(["malgunbd.ttf", "malgun.ttf"], 30)
F_TEXT = _font(["malgun.ttf"], 20)
F_SMALL = _font(["malgun.ttf"], 15)
F_MONO = _font(["consola.ttf", "cour.ttf"], 20)
F_BASE = _font(["consolab.ttf", "consola.ttf"], 18)


def _wrap(draw, text, font, width):
    """글자 폭을 재서 width 안에 들어가도록 줄바꿈"""
    lines = []
    for para in text.split("\n"):
        cur = ""
        for ch in para:
            if draw.textlength(cur + ch, font=font) > width and cur:
                lines.append(cur)
                cur = ch
            else:
                cur += ch
        lines.append(cur)
    return lines


def _dashed_hline(draw, x0, x1, y, color, dash=6):
    for x in range(int(x0), int(x1), dash * 2):
        draw.line([(x, y), (min(x + dash, x1), y)], fill=color, width=1)


class _Canvas:
    """세로로 내용을 이어 붙이는 그리기 도우미 (높이는 마지막에 잘라냄)"""

    def __init__(self, max_h=20000):
        self.img = Image.new("RGB", (W, max_h), "white")
        self.d = ImageDraw.Draw(self.img)
        self.y = M

    def text(self, text, font=F_TEXT, color=(20, 20, 20), gap=6, bg=None):
        for line in _wrap(self.d, text, font, W - 2 * M):
            box = self.d.textbbox((M, self.y), line or " ", font=font)
            h = max(box[3] - self.y, font.size)
            if bg:
                self.d.rectangle([M - 8, self.y - 3, W - M + 8, self.y + h + 5], fill=bg)
            self.d.text((M, self.y), line, font=font, fill=color)
            self.y += h + gap

    def space(self, h):
        self.y += h

    def rule(self):
        self.d.line([(M, self.y), (W - M, self.y)], fill=(210, 210, 210), width=2)
        self.y += 14

    def finish(self):
        return self.img.crop((0, 0, W, self.y + M))


# ---------- 파일 하나 ----------
# 정렬 줄 스타일: (글꼴, 글자색, 바탕색)
F_MONO_B = _font(["consolab.ttf", "consola.ttf"], 20)
F_TEXT_B = _font(["malgunbd.ttf", "malgun.ttf"], 20)
_STYLE = {
    "q": (F_MONO, (23, 78, 166), (232, 240, 254)),      # 찾을 서열: 파랑
    "r": (F_MONO, (19, 115, 51), (230, 244, 234)),      # 결과 서열: 초록
    "mm": (F_MONO_B, (255, 255, 255), (217, 48, 37)),   # 불일치: 빨간 바탕 흰 글자
    "ins": (F_MONO_B, (255, 255, 255), (132, 48, 206)),  # 그 뒤에 삽입: 보라 바탕
    "bar": (F_MONO, (154, 160, 166), None),
    "mark": (F_MONO_B, (217, 48, 37), None),
    "num": (F_MONO, (128, 134, 139), None),
    "head": (F_TEXT_B, (20, 20, 20), None),
    "note": (F_TEXT, (95, 99, 104), None),
    "mmtext": (F_TEXT, (217, 48, 37), None),
    "oktext": (F_TEXT, (19, 115, 51), None),
    "instext": (F_TEXT, (132, 48, 206), None),
    "": (F_TEXT, (20, 20, 20), None),
}
_LABEL_COLOR = {"q": (23, 78, 166), "r": (19, 115, 51)}


def _alignment_block(c, lines):
    """
    alignment_view.build_lines 의 줄 목록을 그림.
    왼쪽 = 서열이 나온 파일 이름, 오른쪽 = 색칠된 서열 (불일치는 빨간 바탕)
    """
    label_w = max([c.d.textlength(l["label"], font=F_TEXT_B) for l in lines if l.get("label")] or [0]) + 24
    x0 = M + label_w
    line_h = 28
    for l in lines:
        if not l["cells"] and not l.get("label"):
            c.y += 10
            continue
        if l.get("label"):
            c.d.text((M, c.y - 2), l["label"], font=F_TEXT_B, fill=_LABEL_COLOR.get(l.get("label_tag"), (60, 60, 60)))
        x = x0
        for text, tag in l["cells"]:
            font, fg, bg = _STYLE.get(tag, _STYLE[""])
            for ch in text:
                w = c.d.textlength(ch, font=font)
                if x + w > W - M:          # 긴 설명은 서열 시작 위치에서 줄바꿈
                    x = x0
                    c.y += line_h
                if bg and ch != " ":
                    c.d.rectangle([x, c.y, x + w, c.y + line_h - 2], fill=bg)
                c.d.text((x, c.y + 2), ch, font=font, fill=fg)
                x += w
        c.y += line_h


def _overview(c, tr, q, r):
    """판독 전체 품질 막대 (깨끗한 구간 / 찾은 위치 표시)"""
    c.text("전체 품질  (회색 막대 = 염기별 품질 Q, 초록 = 깨끗한 구간, 주황 = 찾은 위치, 빨간 점선 = Q20)", F_SMALL,
           (80, 80, 80), gap=4)
    h, x0, w = 90, M, W - 2 * M
    y0 = c.y
    n = len(tr)
    sx = w / n
    c.d.rectangle([x0, y0, x0 + w, y0 + h], outline=(200, 200, 200))
    if q and q["clean"]:
        c.d.rectangle([x0 + (q["clean"][0] - 1) * sx, y0 + 1, x0 + q["clean"][1] * sx, y0 + h - 1],
                      fill=(227, 244, 230))
    if r:
        c.d.rectangle([x0 + (r["start"] - 1) * sx, y0 + 1, x0 + r["end"] * sx, y0 + h - 1],
                      fill=(255, 216, 168), outline=(240, 140, 0))
    for i, qv in enumerate(tr.qual):
        x = x0 + i * sx
        c.d.line([(x, y0 + h - 1), (x, y0 + h - 1 - min(qv, 62) / 62 * (h - 20))], fill=(150, 155, 160))
    _dashed_hline(c.d, x0, x0 + w, y0 + h - 1 - GOOD_Q / 62 * (h - 20), (224, 49, 49))
    for k in range(0, n, 100):
        c.d.text((x0 + k * sx + 2, y0 + 2), str(k or 1), font=F_SMALL, fill=(90, 90, 90))
    c.y += h + 20


def _chromatogram(c, tr, r, compact=False):
    """찾은 위치 주변 크로마토그램 (BASES_PER_ROW 염기씩 여러 줄). compact=True 면 높이를 줄여서 그림"""
    if not compact:
        c.text("크로마토그램  (연한 주황 = 찾은 서열 구간, 노랑 = 차이 있는 염기, 아래 막대 = 품질 Q, 빨간 점선 = Q20)",
               F_SMALL, (80, 80, 80), gap=4)
    s_all = max(1, r["start"] - FLANK)
    e_all = min(len(tr), r["end"] + FLANK)
    diff_pos = set(r.get("diff_pos", []))
    trace_h, base_gap, qbar_h = (110, 26, 36) if compact else (170, 26, 60)
    row_h = 22 + trace_h + base_gap + qbar_h + (16 if compact else 30)
    # 줄 수를 정한 뒤 염기를 줄마다 고르게 나눔 (예: 67개 -> 1줄, 160개 -> 2줄 x 80)
    total = e_all - s_all + 1
    n_rows = -(-total // BASES_PER_ROW)
    per_row = -(-total // n_rows)
    # 가로 비율은 모든 줄이 같게: 평균 피크 간격 x per_row 가 이미지 폭에 맞도록
    gaps = sorted(tr.ploc[i] - tr.ploc[i - 1] for i in range(s_all, e_all))
    spacing = gaps[len(gaps) // 2] if gaps else 12
    px = (W - 2 * M) / (spacing * per_row + 24)
    for s in range(s_all, e_all + 1, per_row):
        e = min(e_all, s + per_row - 1)
        y0 = c.y
        x_lo = tr.ploc[s - 1] - 8
        x_hi = tr.ploc[e - 1] + 8
        # 피크 간격이 평소보다 넓은 줄은 폭을 넘지 않게 그 줄만 조금 줄임
        row_px = min(px, (W - 2 * M) / max(x_hi - x_lo, 1))
        X = lambda scan, px=row_px: M + (scan - x_lo) * px
        top = y0 + 22
        # 배경
        for i in range(s, e + 1):
            xm = X(tr.ploc[i - 1])
            half = 9 * row_px if row_px < 1.5 else 7 * row_px
            if i in diff_pos:
                c.d.rectangle([xm - half, top, xm + half, top + trace_h + base_gap + qbar_h], fill=(255, 224, 102))
            elif r["start"] <= i <= r["end"]:
                c.d.rectangle([xm - half, top, xm + half, top + trace_h + base_gap + qbar_h], fill=(255, 244, 230))
        # 신호
        peak = max(max(ch[x_lo:x_hi + 1] or [1]) for ch in tr.channels.values()) or 1
        for base, ch in tr.channels.items():
            pts = [(X(sc), top + trace_h - ch[sc] / peak * (trace_h - 6))
                   for sc in range(max(0, x_lo), min(len(ch), x_hi + 1))]
            if len(pts) > 1:
                c.d.line(pts, fill=BASE_COLORS.get(base, (0, 0, 0)), width=2, joint="curve")
        # 염기, 번호, 품질
        qtop = top + trace_h + base_gap
        for i in range(s, e + 1):
            xm = X(tr.ploc[i - 1])
            b = tr.seq[i - 1]
            c.d.text((xm, top + trace_h + 4), b, font=F_BASE, fill=BASE_COLORS.get(b, (0, 0, 0)), anchor="mt")
            if i % 10 == 0:
                c.d.text((xm, y0), str(i), font=F_SMALL, fill=(90, 90, 90), anchor="mt")
            qv = tr.qual[i - 1]
            col = (116, 198, 157) if qv >= GOOD_Q else (255, 135, 135)
            c.d.rectangle([xm - 5, qtop + qbar_h - min(qv, 62) / 62 * qbar_h, xm + 5, qtop + qbar_h], fill=col)
        _dashed_hline(c.d, M, X(x_hi), qtop + qbar_h - GOOD_Q / 62 * qbar_h, (224, 49, 49))
        c.d.line([(M, qtop + qbar_h), (X(x_hi), qtop + qbar_h)], fill=(200, 200, 200))
        c.y = y0 + row_h


def render_report(f, r, query, query_name="찾을 서열", frame_setting="자동"):
    """파일 하나의 결과 이미지 (PIL.Image). query_name = 찾을 서열이 나온 파일 이름"""
    c = _Canvas()
    c.text(f"{Path(f).name}", F_TITLE, gap=10)
    c.text(f"찾을 서열: {query_name} ({len(query)} bp)  {query}", F_SMALL, (90, 90, 90))
    c.text(f"분석: {datetime.now():%Y-%m-%d %H:%M}   ·   시퀀싱 서열 검색기", F_SMALL, (140, 140, 140))
    c.space(6)
    c.rule()
    if r["verdict"] == "오류":
        c.text(f"오류: {r['error']}", bg=VERDICT_BG["오류"])
        return c.finish()

    found = r["verdict"] != "없음"
    bg = VERDICT_BG.get(r["verdict"], DIFF_BG)
    if found:
        line = f"판정: {r['verdict']}   |   방향: {r['strand']}   |   위치: {r['start']}-{r['end']}"
        if r.get("site_label"):
            line += f"   |   찾은 부위 품질: {r['site_label']}"
    else:
        line = f"판정: 없음   (가장 비슷한 곳도 {r['diff']}개 차이)"
    c.text(line, F_TEXT, bg=bg, gap=14)

    if found:
        lines, _ = build_lines(r, query, r["read"], query_name, Path(f).name, frame_setting)
        _alignment_block(c, lines)
        if r.get("diff_notes") or r["diffs"]:
            c.text("DNA 차이: " + "; ".join(r.get("diff_notes") or r["diffs"]), F_TEXT, (160, 30, 30))
    q = r.get("qual")
    if q:
        clean = f"{q['clean'][0]}-{q['clean'][1]}" if q["clean"] else "-"
        c.text(f"전체 품질: {q['grade']}   |   깨끗한 구간(Q20 기준): {clean} ({q['clean_len']} bp, 평균 Q{q['clean_meanq']})"
               f"   |   {q['note']}", F_TEXT)
        if found and r["site"]["mixed"]:
            c.text("찾은 부위 이중 피크: " + ", ".join(
                f"{p}번 {b}/{b2} {int(x * 100)}%" for p, b, b2, x in r["site"]["mixed"]), F_TEXT, (160, 30, 30))
    else:
        c.text("(텍스트 파일이라 크로마토그램 품질은 분석하지 않았습니다)", F_SMALL, (120, 120, 120))
    c.space(10)

    tr = r.get("trace")
    if tr is not None:
        _overview(c, tr, q, r if found else None)
        if found:
            _chromatogram(c, tr, r)
    return c.finish()


# ---------- 요약 표 ----------
def _summary_table(c, results):
    """결과 요약 표 (파일마다 한 줄)"""
    cols = [("파일", 260), ("판정", 110), ("위치", 140), ("찾은 부위 품질", 220),
            ("깨끗한 구간", 160), ("전체 품질", 120), ("차이 / 품질 설명", W - 2 * M - 1010)]
    row_pad = 10

    def row(values, bg, font=F_SMALL, bold=False):
        cells = [_wrap(c.d, str(v), font, w - 16) for v, (_, w) in zip(values, cols)]
        h = max(len(lines) for lines in cells) * (font.size + 6) + row_pad * 2
        c.d.rectangle([M, c.y, W - M, c.y + h], fill=bg, outline=(220, 220, 220))
        x = M
        for lines, (_, w) in zip(cells, cols):
            yy = c.y + row_pad
            for ln in lines:
                c.d.text((x + 8, yy), ln, font=F_TEXT if bold else font, fill=(20, 20, 20))
                yy += font.size + 6
            x += w
        c.y += h

    row([n for n, _ in cols], (233, 236, 239), bold=True)
    for f, r in results.items():
        if r["verdict"] == "오류":
            row([Path(f).name, "오류", "", "", "", "", r["error"]], VERDICT_BG["오류"])
            continue
        found = r["verdict"] != "없음"
        q = r.get("qual")
        pos = f"{r['start']}-{r['end']}" + (" (역)" if r["strand"] == "역상보" else "") if found else ""
        clean = (f"{q['clean'][0]}-{q['clean'][1]}" if q["clean"] else "없음") if q else "-"
        if not found:
            detail = f"가장 비슷한 곳도 {r['diff']}개 차이"
        else:
            detail = "; ".join(r.get("diff_notes") or r["diffs"])
        if q:
            detail = (detail + " | " if detail else "") + q["note"]
        row([Path(f).name, r["verdict"], pos, r.get("site_label", "-") if found else "",
             clean, q["grade"] if q else "-", detail],
            VERDICT_BG.get(r["verdict"], DIFF_BG))


def _header(c, title, results, query, query_name):
    c.text(title, F_TITLE, gap=10)
    c.text(f"찾을 서열: {query_name} ({len(query)} bp)  {query}", F_SMALL, (90, 90, 90))
    ok = sum(r["verdict"] == "일치" for r in results.values())
    c.text(f"분석: {datetime.now():%Y-%m-%d %H:%M}   ·   {len(results)}개 중 일치 {ok}개   ·   시퀀싱 서열 검색기",
           F_SMALL, (140, 140, 140))
    c.space(8)


def render_summary(results, query, query_name="찾을 서열"):
    """results: {파일경로: 결과 dict} -> 요약 표 이미지"""
    c = _Canvas()
    _header(c, "시퀀싱 분석 결과 요약", results, query, query_name)
    _summary_table(c, results)
    return c.finish()


def render_overview(results, query, query_name="찾을 서열", frame_setting="자동"):
    """
    모아보기: 요약 표 + 모든 파일의 DNA/아미노산 정렬을 한 표에 + 파일별 크로마토그램
    """
    c = _Canvas(max_h=60000)
    _header(c, "시퀀싱 결과 모아보기", results, query, query_name)
    _summary_table(c, results)
    c.space(24)
    c.rule()
    _alignment_block(c, combined_dna_lines(results, query, query_name))
    c.space(10)
    _alignment_block(c, combined_aa_lines(results, query, query_name, frame_setting))
    with_trace = [(f, r) for f, r in results.items()
                  if r.get("trace") is not None and r["verdict"] not in ("없음", "오류")]
    if with_trace:
        c.space(10)
        c.rule()
        c.text("크로마토그램 모아보기  (연한 주황 = 찾은 서열 구간, 노랑 = 차이 있는 염기, 아래 막대 = 품질 Q, 빨간 점선 = Q20)",
               F_TEXT_B, gap=10)
        for f, r in with_trace:
            c.text(f"{Path(f).name}   ·   {r['verdict']}   ·   위치 {r['start']}-{r['end']}"
                   + (f"   ·   찾은 부위 품질 {r['site_label']}" if r.get("site_label") else ""),
                   F_TEXT, bg=VERDICT_BG.get(r["verdict"], DIFF_BG), gap=4)
            _chromatogram(c, r["trace"], r, compact=True)
            c.space(6)
    return c.finish()


def export_pdf(path, results, query, query_name="찾을 서열", frame_setting="자동"):
    """PDF 한 파일: 1쪽 = 모아보기, 이후 = 파일별 상세"""
    pages = [render_overview(results, query, query_name, frame_setting)]
    pages += [render_report(f, r, query, query_name, frame_setting) for f, r in results.items()]
    pages[0].save(path, "PDF", resolution=150, save_all=True, append_images=pages[1:])
