"""
시퀀싱 결과(R115A-1_P1.txt 등)에 기준 서열(Primer_..._R115A_F.txt 등)이
들어 있는지 확인하는 프로그램.

사용법:
    python analyze_seq.py            # 이 파일이 있는 폴더 분석
    python analyze_seq.py <폴더경로>  # 다른 폴더 분석

결과는 화면에 출력되고, 같은 폴더에 analysis_result.csv 로도 저장됩니다.
"""
import csv
import re
import sys
from pathlib import Path

# 시퀀싱 결과 파일 이름 규칙: <변이이름>-<반복번호>_<시퀀싱프라이머>.txt  예) R115A-1_P1.txt
SAMPLE_PATTERN = re.compile(r"^(?P<mut>[A-Za-z]\d+[A-Za-z])-(?P<rep>\d+)_(?P<sp>[^.]+)\.txt$")

# 차이(치환+삽입+결실) 비율이 이 값 이하이면 '서열이 있음(차이 있음)'으로 판단
MAX_DIFF_RATIO = 0.2


def clean_seq(text):
    """FASTA 헤더(>로 시작하는 줄)를 빼고 서열만 이어 붙여 대문자로 반환"""
    lines = [l.strip() for l in text.splitlines() if l.strip() and not l.strip().startswith(">")]
    return re.sub(r"[^A-Za-z]", "", "".join(lines)).upper()


def read_seq(path):
    return clean_seq(Path(path).read_text(encoding="utf-8", errors="ignore"))


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def align(query, read):
    """
    query(기준 서열) 전체가 read(시퀀싱 결과)의 어디에 가장 잘 맞는지 찾는다.
    치환/삽입/결실을 모두 1개의 차이로 계산 (semi-global 정렬).
    반환: (차이 개수, 시작 위치(1부터), 끝 위치, 기준 정렬문자열, 표시줄, 결과 정렬문자열)
    """
    m, n = len(query), len(read)
    # dp[i][j]: query 앞 i글자를 read 의 j번째에서 끝나도록 맞출 때 최소 차이
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        dp[i][0] = i
        qi = query[i - 1]
        prev, cur = dp[i - 1], dp[i]
        for j in range(1, n + 1):
            match = 0 if (qi == read[j - 1] and qi != "N") else 1
            cur[j] = min(prev[j - 1] + match, prev[j] + 1, cur[j - 1] + 1)

    end = min(range(n + 1), key=lambda j: dp[m][j])
    diff = dp[m][end]

    # 역추적해서 정렬 모양 만들기
    i, j = m, end
    q_aln, r_aln = [], []
    while i > 0:
        if j > 0 and dp[i][j] == dp[i - 1][j - 1] + (0 if query[i - 1] == read[j - 1] != "N" else 1):
            q_aln.append(query[i - 1]); r_aln.append(read[j - 1]); i -= 1; j -= 1
        elif dp[i][j] == dp[i - 1][j] + 1:
            q_aln.append(query[i - 1]); r_aln.append("-"); i -= 1   # 결과에서 빠짐(결실)
        else:
            q_aln.append("-"); r_aln.append(read[j - 1]); j -= 1    # 결과에 추가됨(삽입)
    q_aln, r_aln = "".join(reversed(q_aln)), "".join(reversed(r_aln))
    bar = "".join("|" if a == b else " " for a, b in zip(q_aln, r_aln))
    return diff, j + 1, end, q_aln, bar, r_aln


def describe_diffs(q_aln, r_aln):
    """차이를 '12:A>G', '20:+T', '31:C>-' 형식(기준 서열 위치 기준)으로 정리"""
    out, pos = [], 0
    for a, b in zip(q_aln, r_aln):
        if a != "-":
            pos += 1
        if a == b:
            continue
        if a == "-":
            out.append(f"{pos}:+{b}")       # pos 다음에 삽입
        elif b == "-":
            out.append(f"{pos}:{a}>-")      # 결실
        else:
            out.append(f"{pos}:{a}>{b}")    # 치환
    return out


def diff_read_positions(q_aln, r_aln, start):
    """각 차이가 결과 서열의 몇 번째 염기(1부터)에 해당하는지 (describe_diffs 와 같은 순서)"""
    out, rpos = [], start - 1
    for a, b in zip(q_aln, r_aln):
        if b != "-":
            rpos += 1
        if a != b:
            # 결실(-)이면 결과 서열에 염기가 없으므로 바로 앞 염기 위치
            out.append(max(rpos, 1))
    return out


def analyze_one(ref, read, max_diff_ratio=MAX_DIFF_RATIO):
    """ref 가 read 안에 있는지 정방향/역상보 모두 검사해서 결과 dict 반환"""
    fwd = align(ref, read)
    rev = align(revcomp(ref), read)
    strand, res = ("정방향", fwd) if fwd[0] <= rev[0] else ("역상보", rev)
    diff, start, end, q_aln, bar, r_aln = res
    if diff == 0:
        verdict = "일치"
    elif diff <= len(ref) * max_diff_ratio:
        verdict = f"차이 {diff}개"
    else:
        verdict = "없음"
    return {
        "verdict": verdict, "diff": diff, "strand": strand, "start": start, "end": end,
        "q_aln": q_aln, "bar": bar, "r_aln": r_aln, "diffs": describe_diffs(q_aln, r_aln),
        "diff_pos": diff_read_positions(q_aln, r_aln, start),
    }


def find_reference(folder, mut, sample_files):
    """변이 이름(예: R115A)이 들어간, 시퀀싱 결과가 아닌 txt 파일을 기준 서열로 사용"""
    cands = [p for p in folder.glob("*.txt")
             if re.search(rf"(^|[_\-]){re.escape(mut)}([_\-.]|$)", p.name)
             and p not in sample_files]
    return sorted(cands)


def main():
    folder = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent

    samples = {}
    for p in sorted(folder.glob("*.txt")):
        m = SAMPLE_PATTERN.match(p.name)
        if m:
            samples.setdefault(m["mut"], []).append(p)
    all_sample_files = {p for ps in samples.values() for p in ps}

    if not samples:
        print("시퀀싱 결과 파일(예: R115A-1_P1.txt)을 찾지 못했습니다.")
        return

    rows = []
    for mut, files in samples.items():
        refs = find_reference(folder, mut, all_sample_files)
        print("=" * 70)
        if not refs:
            print(f"[{mut}] 기준 서열 파일을 찾지 못했습니다. (이름에 '{mut}'가 들어간 파일 필요)")
            continue
        if len(refs) > 1:
            print(f"[{mut}] 기준 서열 후보가 여러 개입니다: {[r.name for r in refs]} -> 첫 번째 사용")
        ref_path = refs[0]
        ref = read_seq(ref_path)
        print(f"[{mut}] 기준 서열: {ref_path.name} ({len(ref)} bp)")
        print(f"        {ref}")

        for f in files:
            read = read_seq(f)
            r = analyze_one(ref, read)
            verdict, diff, strand, start, end = r["verdict"], r["diff"], r["strand"], r["start"], r["end"]
            q_aln, bar, r_aln, diffs = r["q_aln"], r["bar"], r["r_aln"], r["diffs"]

            print("-" * 70)
            print(f"  {f.name}  (서열 길이 {len(read)} bp)")
            print(f"  판정: {verdict}  |  방향: {strand}  |  위치: {start}-{end}")
            if verdict != "없음":
                print(f"    기준 {q_aln}")
                print(f"         {bar}")
                print(f"    결과 {r_aln}")
                if diffs:
                    print(f"    차이 위치(기준 서열 기준): {', '.join(diffs)}")

            rows.append({
                "변이": mut, "파일": f.name, "기준서열파일": ref_path.name,
                "판정": verdict, "차이개수": diff, "방향": strand,
                "시작": start, "끝": end, "차이상세": "; ".join(diffs),
            })

    out = folder / "analysis_result.csv"
    with open(out, "w", newline="", encoding="utf-8-sig") as fp:   # 엑셀에서 한글 안 깨지게
        w = csv.DictWriter(fp, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("=" * 70)
    print(f"결과 저장: {out}")


if __name__ == "__main__":
    main()
