# SeqFinder – 시퀀싱 서열 검색기

Sanger 시퀀싱 결과 파일에서 원하는 서열(예: 돌연변이 프라이머 서열)이 들어 있는지 확인하고,
`.ab1` 원본 파일의 크로마토그램 품질까지 함께 분석하는 Windows 프로그램입니다.

## 주요 기능

- **서열 검색**: 찾을 서열을 입력하면 시퀀싱 결과에서 위치를 찾습니다.
  - 정방향 / 역상보(reverse complement) 모두 검색
  - 치환·삽입·결실을 허용하는 정렬로, 차이가 있으면 정확한 위치 표시
- **크로마토그램 품질 분석** (`.ab1` 파일)
  - 찾은 부위의 품질 점수(Q)와 이중 피크 여부
  - 깨끗한 구간(Q20 기준) 범위와 길이
  - 가운데에서 품질이 급락하거나 두 서열이 겹쳐 읽히는 경우 감지
  - 서열 차이가 판독 오류인지 판단할 근거 표시 (품질, 피크 간격, 이중 피크)
- **결과 저장**: CSV(엑셀), PNG/JPG 이미지(파일별 상세 + 전체 요약 표)

## 사용 방법 (exe)

1. [Releases](../../releases)에서 `SeqFinder.exe` 를 내려받아 실행합니다. (Python 설치 불필요, Windows 10/11)
2. **찾을 서열** 칸에 서열을 붙여넣거나 `파일에서 불러오기` 를 누릅니다.
3. **파일 추가** 로 시퀀싱 결과 파일을 여러 개 선택합니다. (`.ab1`, `.txt`, `.seq`, `.fasta`)
4. **▶ 분석 시작** → 결과 표의 행을 클릭하면 정렬과 크로마토그램이 표시됩니다.
5. `CSV로 저장`, `선택 결과 이미지 저장`, `전체 이미지 저장` 으로 결과를 저장합니다.

> 서명되지 않은 exe라서 Windows SmartScreen 경고가 나오면 **추가 정보 → 실행** 을 누르세요.

## 판정 기준

| 항목 | 기준 | 설정 위치 |
|---|---|---|
| 서열 판정 | 차이 0개 = 일치 / 찾을 서열 길이의 20% 이하 = 차이 N개 / 초과 = 없음 | 화면의 `허용 차이 비율(%)`, `analyze_seq.py` `MAX_DIFF_RATIO` |
| 믿을 수 있는 염기 | 품질 Q20 이상 (오류율 1%) | `abi_reader.py` `GOOD_Q` |
| 이중 피크 의심 | 두 번째 피크가 첫 번째 피크의 35% 이상 | `abi_reader.py` `MIXED_RATIO` |
| 깨끗한 구간 | Mott 트리밍 (Phred/Biopython 과 같은 방식, Q20 기준) | `abi_reader.py` `clean_region()` |
| 전체 품질 "좋음" | 깨끗한 구간 700 bp 이상 | `abi_reader.py` `GOOD_LEN` |

시퀀싱 결과의 앞쪽·뒤쪽 몇십 염기가 불안정한 것은 정상이며, 판정은 **가운데 깨끗한 구간**을 기준으로 합니다.

## 소스에서 실행하기

Python 3.11 에서 테스트했습니다.

```bash
pip install -r requirements.txt
python seq_finder_gui.py
```

명령줄 버전(폴더 안의 `R115A-1_P1.txt` 같은 결과 파일을, 이름에 같은 변이 이름이 들어간 기준 서열 파일과 자동으로 짝지어 비교):

```bash
python analyze_seq.py <폴더경로>
```

## exe 다시 만들기

```bash
pip install pyinstaller
python -m PyInstaller --noconfirm --onefile --windowed --name SeqFinder --exclude-module numpy seq_finder_gui.py
```

결과물은 `dist/SeqFinder.exe` 에 생깁니다.
`--exclude-module numpy` 를 빼면 PC에 설치된 numpy 가 딸려 들어가 exe 가 수백 MB 로 커질 수 있습니다.

## 파일 구성

| 파일 | 내용 |
|---|---|
| `seq_finder_gui.py` | 프로그램 화면 (tkinter) |
| `analyze_seq.py` | 서열 정렬·비교 로직, 명령줄 버전 |
| `abi_reader.py` | `.ab1`(ABIF) 파일 읽기와 크로마토그램 품질 분석 (외부 라이브러리 없음) |
| `report_image.py` | 결과 이미지 그리기 (Pillow) |

## 데이터 파일에 대해

`.gitignore` 에서 시퀀싱 데이터(`.ab1`, `.pdf`, 결과 `.txt`, 프라이머 서열 파일)는 저장소에 올라가지 않도록 막아 두었습니다.
미발표 연구 데이터가 실수로 공개되지 않게 하기 위해서입니다.
