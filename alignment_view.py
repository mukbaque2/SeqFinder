"""
정렬 결과를 '보여주기 좋은 줄' 목록으로 만드는 모듈.
화면(seq_finder_gui.py)과 이미지(report_image.py)가 같은 내용을 그리도록 공통으로 사용합니다.

각 줄은 {"label": 왼쪽 이름표, "label_tag": 이름표 스타일, "cells": [(글자, 스타일), ...]} 형태입니다.
스타일 이름:
  head   : 소제목            note : 설명 문장
  q      : 찾을 서열 글자      r    : 결과 서열 글자
  mm     : 일치하지 않는 글자 (두 줄 모두)
  bar    : 일치 표시 '|'      mark : 불일치 표시 '*'
  num    : 위치 번호
"""
from analyze_seq import revcomp

BLOCK_NT = 60    # DNA 한 줄에 보여줄 글자 수
BLOCK_AA = 60    # 아미노산 한 줄에 보여줄 글자 수

_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AA[16 * i + 4 * j + k]
         for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}
AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
       "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
       "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Stop", "X": "?"}


def translate(seq, frame=1):
    """frame(1~3)부터 3개씩 끊어 번역. N 등이 섞인 코돈은 X"""
    s = seq[frame - 1:]
    return "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def auto_frame(query):
    """
    찾을 서열에서 정지 코돈(*)이 가장 적은 읽는 틀.
    반환: (틀 번호, 정지 코돈이 없는 틀 목록)
    """
    stops = {f: translate(query, f).count("*") for f in (1, 2, 3)}
    best = min(stops, key=lambda f: (stops[f], f))
    return best, [f for f in (1, 2, 3) if stops[f] == 0]


def _seq_cells(seq, other, kind):
    """seq 의 글자마다 스타일 지정: other 와 다르면 mm"""
    return [(a, "mm" if a != b else kind) for a, b in zip(seq, other)]


def _num(n, width):
    return (str(n) if n is not None else "").rjust(width) + " "


def dna_lines(r, query_name, read_name):
    """DNA 정렬 줄 목록"""
    q, rd = r["q_aln"], r["r_aln"]
    rev = r["strand"] == "역상보"
    q_len = sum(c != "-" for c in q)
    width = len(str(max(q_len, r["end"])))
    lines = [{"label": "", "cells": [("DNA 정렬", "head"),
                                     (f"   (위치 {r['start']}-{r['end']}, {r['strand']})", "note")]}]
    # 찾을 서열 번호: 역상보면 뒤에서부터 센다
    q_done, rpos = 0, r["start"] - 1
    q_label = query_name + (" (역상보)" if rev else "")
    for i in range(0, len(q), BLOCK_NT):
        qs, rs = q[i:i + BLOCK_NT], rd[i:i + BLOCK_NT]
        nq = sum(c != "-" for c in qs)
        nr = sum(c != "-" for c in rs)
        if nq:
            q_from, q_to = ((q_len - q_done, q_len - q_done - nq + 1) if rev
                            else (q_done + 1, q_done + nq))
        else:
            q_from = q_to = None
        q_done += nq
        r_from, rpos = (rpos + 1 if nr else None), rpos + nr
        r_to = rpos if nr else None
        marks = [("|", "bar") if a == b else ("*", "mark") for a, b in zip(qs, rs)]
        lines.append({"label": q_label, "label_tag": "q",
                      "cells": [(_num(q_from, width), "num")] + _seq_cells(qs, rs, "q")
                      + [(" " + _num(q_to, 0), "num")]})
        lines.append({"label": "", "cells": [(" " * (width + 1), "num")] + marks})
        lines.append({"label": read_name, "label_tag": "r",
                      "cells": [(_num(r_from, width), "num")] + _seq_cells(rs, qs, "r")
                      + [(" " + _num(r_to, 0), "num")]})
        lines.append({"label": "", "cells": []})
    return lines


def aa_compare(query, read, r, frame):
    """
    찾을 서열과 결과 서열(찾은 구간)을 같은 읽는 틀로 번역해서 비교.
    결과 서열은 찾을 서열과 같은 방향으로 맞춘 뒤 번역합니다.
    """
    seg = read[r["start"] - 1:r["end"]]
    if r["strand"] == "역상보":
        seg = revcomp(seg)
    qa, ra = translate(query, frame), translate(seg, frame)
    diffs = []
    for k in range(max(len(qa), len(ra))):
        a = qa[k] if k < len(qa) else "-"
        b = ra[k] if k < len(ra) else "-"
        if a != b:
            qc = query[frame - 1 + 3 * k: frame + 2 + 3 * k]
            rc = seg[frame - 1 + 3 * k: frame + 2 + 3 * k]
            diffs.append((k + 1, a, b, qc, rc))
    shift = (len(seg) - len(query)) % 3 != 0
    return {"frame": frame, "q": qa, "r": ra, "diffs": diffs, "frameshift": shift}


def aa_lines(aa, query_name, read_name, frame_note=""):
    """아미노산 정렬 줄 목록"""
    qa, ra = aa["q"], aa["r"]
    n = max(len(qa), len(ra))
    qa, ra = qa.ljust(n, "-"), ra.ljust(n, "-")
    width = len(str(n))
    head = f"   (읽는 틀 +{aa['frame']}{frame_note})"
    lines = [{"label": "", "cells": [("아미노산 정렬", "head"), (head, "note")]}]
    for i in range(0, n, BLOCK_AA):
        qs, rs = qa[i:i + BLOCK_AA], ra[i:i + BLOCK_AA]
        end = i + len(qs)
        marks = [("|", "bar") if a == b else ("*", "mark") for a, b in zip(qs, rs)]
        lines.append({"label": query_name, "label_tag": "q",
                      "cells": [(_num(i + 1, width), "num")] + _seq_cells(qs, rs, "q") + [(" " + str(end), "num")]})
        lines.append({"label": "", "cells": [(" " * (width + 1), "num")] + marks})
        lines.append({"label": read_name, "label_tag": "r",
                      "cells": [(_num(i + 1, width), "num")] + _seq_cells(rs, qs, "r") + [(" " + str(end), "num")]})
    if aa["diffs"]:
        desc = [f"{k}번 {AA3.get(a, a)}({a})→{AA3.get(b, b)}({b}) [{qc or '-'}→{rc or '-'}]"
                for k, a, b, qc, rc in aa["diffs"][:12]]
        more = f" 외 {len(aa['diffs']) - 12}개" if len(aa["diffs"]) > 12 else ""
        lines.append({"label": "", "cells": [("아미노산 차이: " + ", ".join(desc) + more, "mmtext")]})
        if aa["frameshift"]:
            lines.append({"label": "", "cells": [
                ("※ 염기 수가 3의 배수만큼 차이 나지 않아 이후 아미노산이 모두 바뀝니다 (프레임시프트). "
                 "판독 오류인지 크로마토그램을 확인하세요.", "mmtext")]})
    else:
        lines.append({"label": "", "cells": [("아미노산 서열이 모두 일치합니다.", "oktext")]})
    lines.append({"label": "", "cells": []})
    return lines


def frame_note(query, setting):
    """
    읽는 틀 결정. setting: '자동' 또는 '1'/'2'/'3'
    반환: (틀 번호, 화면에 붙일 설명)
    """
    if setting in ("1", "2", "3"):
        return int(setting), ", 직접 선택"
    f, clean = auto_frame(query)
    if len(clean) == 1:
        return f, ", 자동: 정지 코돈이 없는 유일한 틀"
    if len(clean) > 1:
        return f, (f", 자동: 정지 코돈 없는 틀이 {'/'.join(map(str, clean))}번 여러 개 → "
                   "맞는 틀인지 확인하고 필요하면 '읽는 틀'에서 직접 고르세요")
    return f, ", 자동: 모든 틀에 정지 코돈이 있음 → 확인 필요"


def build_lines(r, query, read, query_name, read_name, frame_setting="자동"):
    """DNA + 아미노산 정렬 전체 줄 목록과 아미노산 비교 결과"""
    lines = dna_lines(r, query_name, read_name)
    frame, note = frame_note(query, frame_setting)
    aa = aa_compare(query, read, r, frame)
    lines += aa_lines(aa, query_name, read_name, note)
    return lines, aa


# ---------- 여러 파일 모아보기 ----------
def _revcomp_aln(s):
    """'-' 가 섞인 정렬 문자열의 역상보 ('-' 는 그대로)"""
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def project_to_query(r):
    """
    결과 서열을 찾을 서열 좌표(1..len)에 맞춰 한 글자씩 대응시킴.
    반환: (글자 목록, {찾을 서열 위치: [그 뒤에 삽입된 글자들]})
    - 결실이면 '-' / 삽입은 바로 앞 위치에 기록 (모아보기 표에서 보라색으로 표시)
    - 역상보로 찾은 경우 찾을 서열 방향으로 뒤집어서 맞춤
    """
    q, rd = r["q_aln"], r["r_aln"]
    if r["strand"] == "역상보":
        q, rd = _revcomp_aln(q), _revcomp_aln(rd)
    cells, ins, k = [], {}, 0
    for a, b in zip(q, rd):
        if a != "-":
            cells.append(b)
            k += 1
        else:
            ins.setdefault(k, []).append(b)
    return cells, ins


def _ruler(n, width, start=1):
    """10 단위 위치 눈금 줄"""
    s = [" "] * n
    for p in range(start, start + n):
        if p % 10 == 0:
            label = str(p)
            i = p - start - len(label) + 1
            for j, ch in enumerate(label):
                if 0 <= i + j < n:
                    s[i + j] = ch
    return [(" " * (width + 1), "num"), ("".join(s), "num")]


def combined_dna_lines(results, query, query_name, block=BLOCK_NT):
    """찾을 서열 1줄 + 결과 파일마다 1줄씩, 같은 위치끼리 세로로 맞춘 DNA 정렬"""
    found = [(f, r) for f, r in results.items() if r.get("verdict") not in ("없음", "오류")]
    proj = {f: project_to_query(r) for f, r in found}
    width = max([len(str(len(query)))] + [len(str(r["end"])) for _, r in found])
    lines = [{"label": "", "cells": [("DNA 정렬 모아보기", "head"),
                                     ("   (찾을 서열 위치 기준, 빨강 = 불일치, 보라 = 그 뒤에 염기 삽입)", "note")]}]
    for i in range(0, len(query), block):
        qs = query[i:i + block]
        lines.append({"label": "", "cells": _ruler(len(qs), width, i + 1)})
        lines.append({"label": query_name, "label_tag": "q",
                      "cells": [(_num(i + 1, width), "num"), (qs, "q"), (" " + str(i + len(qs)), "num")]})
        for f, r in found:
            cells, ins = proj[f]
            row = [(_num(None, width), "num")]
            for k, ch in enumerate(cells[i:i + block]):
                pos = i + k + 1
                if pos in ins:
                    row.append((ch, "ins"))
                else:
                    row.append((ch, "r" if ch == query[pos - 1] else "mm"))
            tail = []
            if i + block >= len(query):   # 마지막 줄에 결과 서열 위치와 삽입 정보
                tail.append((f" {r['start']}-{r['end']}", "num"))
                if r["strand"] == "역상보":
                    tail.append((" (역상보)", "note"))
                extra = [f"{p}번 뒤 +{''.join(b)}" for p, b in sorted(ins.items())]
                if extra:
                    tail.append(("   삽입: " + ", ".join(extra), "instext"))
            lines.append({"label": _fname(f), "label_tag": "r", "cells": row + tail})
        lines.append({"label": "", "cells": []})
    for f, r in results.items():
        if r.get("verdict") in ("없음", "오류"):
            why = r.get("error") or f"찾지 못함 (가장 비슷한 곳도 {r.get('diff')}개 차이)"
            lines.append({"label": _fname(f), "label_tag": "r", "cells": [(why, "mmtext")]})
    return lines


def combined_aa_lines(results, query, query_name, frame_setting="자동", block=BLOCK_AA):
    """찾을 서열 번역 1줄 + 결과 파일마다 번역 1줄"""
    frame, note = frame_note(query, frame_setting)
    found = [(f, r) for f, r in results.items() if r.get("verdict") not in ("없음", "오류")]
    aas = {f: aa_compare(query, r["read"], r, frame) for f, r in found}
    qa = translate(query, frame)
    n = max([len(qa)] + [len(a["r"]) for a in aas.values()])
    width = len(str(n))
    lines = [{"label": "", "cells": [("아미노산 정렬 모아보기", "head"), (f"   (읽는 틀 +{frame}{note})", "note")]}]
    for i in range(0, n, block):
        qs = qa[i:i + block]
        lines.append({"label": "", "cells": _ruler(min(block, n - i), width, i + 1)})
        lines.append({"label": query_name, "label_tag": "q",
                      "cells": [(_num(i + 1, width), "num"), (qs, "q"), (" " + str(i + len(qs)), "num")]})
        for f, _ in found:
            a = aas[f]
            rs = a["r"][i:i + block]
            row = [(_num(None, width), "num")]
            row += [(ch, "r" if k + i < len(qa) and ch == qa[k + i] else "mm") for k, ch in enumerate(rs)]
            tail = ""
            if i + block >= n:
                tail = "  일치" if not a["diffs"] else f"  차이 {len(a['diffs'])}개" + (
                    " (프레임시프트)" if a["frameshift"] else "")
            lines.append({"label": _fname(f), "label_tag": "r",
                          "cells": row + [(tail, "oktext" if not a["diffs"] else "mmtext")]})
        lines.append({"label": "", "cells": []})
    return lines


def _fname(f):
    """경로 -> 파일 이름"""
    return str(f).replace("\\", "/").rsplit("/", 1)[-1]
