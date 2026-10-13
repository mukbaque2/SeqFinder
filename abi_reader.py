"""
ABI 시퀀싱 원본 파일(.ab1, ABIF 형식) 읽기 + 크로마토그램 품질 분석.
외부 라이브러리 없이 동작합니다.
"""
import struct
from pathlib import Path

# ABIF 데이터 형식 코드 -> struct 형식
_TYPES = {1: "B", 2: "s", 3: "H", 4: "h", 5: "i", 7: "f", 8: "d", 18: "p", 19: "s"}

GOOD_Q = 20          # 이 점수 이상이면 믿을 수 있는 염기 (오류율 1%)
MIXED_RATIO = 0.35   # 두 번째 피크가 첫 번째의 35% 이상이면 이중 피크 의심
GOOD_LEN = 700       # 깨끗한 구간이 이 길이 이상이면 정상적인 판독


def read_abif_tags(path):
    """ab1 파일의 모든 태그를 {'PBAS2': 값, ...} 형태로 반환"""
    data = Path(path).read_bytes()
    if data[:4] != b"ABIF":
        raise ValueError("ab1(ABIF) 파일이 아닙니다.")
    # 헤더의 루트 디렉터리 항목: 원소 개수와 위치
    _, _, _, _, n_entries, _, dir_offset, _ = struct.unpack(">4sihhiiii", data[6:34])
    tags = {}
    for k in range(n_entries):
        off = dir_offset + k * 28
        name, number, etype, _, count, size, doff, _ = struct.unpack(">4sihhiiii", data[off:off + 28])
        raw = data[off + 20:off + 24] if size <= 4 else data[doff:doff + size]
        key = name.decode("latin-1") + str(number)
        fmt = _TYPES.get(etype)
        try:
            if etype in (2, 19):
                val = raw[:count].rstrip(b"\x00").decode("latin-1")
            elif etype == 18:
                val = raw[1:1 + raw[0]].decode("latin-1")
            elif fmt:
                val = struct.unpack(f">{count}{fmt}", raw[:count * struct.calcsize(fmt)])
                if count == 1:
                    val = val[0]
            else:
                val = raw
        except (struct.error, IndexError):
            val = raw
        tags[key] = val
    return tags


class Trace:
    """ab1 한 개의 염기서열, 품질 점수, 피크 위치, 4색 신호"""

    def __init__(self, path):
        t = read_abif_tags(path)
        self.path = Path(path)
        self.seq = (t.get("PBAS2") or t.get("PBAS1")).upper()
        q = t.get("PCON2") or t.get("PCON1")
        self.qual = [ord(c) for c in q] if isinstance(q, str) else list(q)
        self.ploc = list(t.get("PLOC2") or t.get("PLOC1"))
        order = t.get("FWO_1", "GATC")
        # DATA9~12 = 분석된 신호, FWO_1 에 적힌 염기 순서대로
        self.channels = {base: list(t[f"DATA{9 + i}"]) for i, base in enumerate(order)}
        self._ratios = None

    def __len__(self):
        return len(self.seq)

    # ---------- 염기 단위 ----------
    def peak_heights(self, i):
        """i번째(0부터) 염기 피크 위치에서 4색 신호 값 {염기: 높이}"""
        x = self.ploc[i]
        return {b: max(ch[max(0, x - 1):x + 2]) for b, ch in self.channels.items()}

    def second_peak(self, i):
        """(두 번째 피크 / 가장 큰 피크 비율, 판독된 염기 외에 가장 큰 염기)"""
        h = self.peak_heights(i)
        top = max(h.values())
        others = {b: v for b, v in h.items() if b != self.seq[i]}
        b2 = max(others, key=others.get)
        vals = sorted(h.values(), reverse=True)
        return (vals[1] / top if top > 0 else 0.0), b2

    @property
    def ratios(self):
        if self._ratios is None:
            self._ratios = [self.second_peak(i)[0] for i in range(len(self.seq))]
        return self._ratios

    def spacing(self, i):
        """i번째 피크와 앞 피크 사이 간격"""
        return self.ploc[i] - self.ploc[i - 1] if i > 0 else None

    # ---------- 구간 단위 ----------
    def clean_region(self, cutoff=0.01):
        """
        Mott 트리밍(Phred/Biopython 과 같은 방식): 오류율 cutoff(=Q20) 기준으로
        품질이 좋은 가장 긴 구간. 반환: (시작, 끝) 1부터 세는 번호, 없으면 None
        """
        best, best_s, best_e = 0.0, None, None
        run, run_s = 0.0, 0
        for i, q in enumerate(self.qual):
            run += cutoff - 10 ** (-q / 10)
            if run <= 0:
                run, run_s = 0.0, i + 1
            elif run > best:
                best, best_s, best_e = run, run_s, i
        return (best_s + 1, best_e + 1) if best_s is not None else None

    def region_stats(self, start, end):
        """[start, end](1부터) 구간의 품질 통계"""
        s, e = max(1, start), min(len(self.seq), end)
        qs = self.qual[s - 1:e]
        low = [(i, self.seq[i - 1], self.qual[i - 1]) for i in range(s, e + 1) if self.qual[i - 1] < GOOD_Q]
        mixed = []
        for i in range(s, e + 1):
            r, b2 = self.second_peak(i - 1)
            if r >= MIXED_RATIO:
                mixed.append((i, self.seq[i - 1], b2, round(r, 2)))
        return {
            "start": s, "end": e, "min_q": min(qs) if qs else 0,
            "mean_q": round(sum(qs) / len(qs), 1) if qs else 0,
            "low": low, "mixed": mixed,
        }

    def _mean(self, arr, s, e):
        part = arr[max(0, s):max(0, e)]
        return sum(part) / len(part) if part else 0

    def summary(self, good_len=GOOD_LEN, window=25):
        """
        파일 전체 품질 요약.
        grade: '좋음'(깨끗한 구간 good_len 이상) / '짧음' / '나쁨'(200bp 미만)
        note : 깨끗한 구간이 끝나는 모양 설명
        """
        region = self.clean_region()
        res = {"length": len(self.seq), "clean": region, "clean_len": 0, "clean_meanq": 0,
               "grade": "나쁨", "note": "깨끗한 구간 없음"}
        if not region:
            return res
        s, e = region
        st = self.region_stats(s, e)
        res.update(clean_len=e - s + 1, clean_meanq=st["mean_q"], clean_mixed=st["mixed"])

        before = self._mean(self.qual, e - window, e)       # 끝나기 직전 품질
        after = self._mean(self.qual, e, e + window)        # 끝난 직후 품질
        mixed_after = self._mean(self.ratios, e, e + 75)    # 끝난 뒤 두 번째 피크 비율
        if res["clean_len"] >= good_len:
            res["grade"], res["note"] = "좋음", f"{e}번 이후 정상적인 끝부분 감소"
        else:
            res["grade"] = "짧음" if res["clean_len"] >= 200 else "나쁨"
            if before >= 40 and after < 20:
                res["note"] = f"{e}번에서 품질 급락"
                if mixed_after >= 0.25:
                    res["note"] += f", 이후 이중 피크(두 서열 겹침, 평균 {int(mixed_after * 100)}%)"
            else:
                res["note"] = f"{e}번 이후 품질 저하"
        return res


def weakest_in_run(tr, i):
    """
    i번째(1부터) 염기가 같은 염기 연속(예: TTT) 안에 있으면 그중 품질이 가장 낮은 위치를 반환.
    (TTT 중 어느 T가 '추가된' 것인지 정렬만으로는 알 수 없어서)
    """
    k = i - 1
    if not (0 <= k < len(tr)):
        return i
    b, s, e = tr.seq[k], k, k
    while s > 0 and tr.seq[s - 1] == b:
        s -= 1
    while e < len(tr) - 1 and tr.seq[e + 1] == b:
        e += 1
    return min(range(s, e + 1), key=lambda j: tr.qual[j]) + 1


def suspicious_reason(tr, i):
    """i번째(1부터) 염기의 품질 정보와 판독 오류 근거 문자열"""
    k = i - 1
    notes = [f"Q{tr.qual[k]}"]
    sp = tr.spacing(k)
    if sp is not None:
        med = sorted(tr.ploc[j + 1] - tr.ploc[j] for j in range(max(0, k - 30), min(len(tr) - 1, k + 30)))
        med = med[len(med) // 2]
        if sp < med * 0.8:
            notes.append(f"피크 간격 좁음({sp}/{med})")
    near = sorted(tr.qual[max(0, k - 10):k] + tr.qual[k + 1:k + 11])
    if near and tr.qual[k] < near[len(near) // 2] * 0.5:
        notes.append(f"주변(Q{near[len(near) // 2]})보다 낮음")
    r, b2 = tr.second_peak(k)
    if r >= MIXED_RATIO:
        notes.append(f"이중 피크 {b2} {int(r * 100)}%")
    if tr.qual[k] < GOOD_Q or len(notes) > 1:
        notes.append("판독 오류 가능")
    return ", ".join(notes)
