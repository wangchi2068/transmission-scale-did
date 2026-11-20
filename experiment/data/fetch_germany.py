"""Reproducible data download script for the German Varicella panel (2014--2019).

Fetches the three public raw datasets into experiment/data/raw/germany/ with LF newlines:
  1. rki_varicella_2014_2019.json
     Source: Robert Koch-Institut (RKI) SurvStat@RKI 2.0 SOAP 1.2 WebService
     Endpoint: https://tools.rki.de/SurvStat/SurvStatWebService.svc
     Query parameters:
       - Cube: SurvStat
       - Pathogen filter: [PathogenOut].[KategorieNz].[Krankheit DE].&[Windpocken]
       - RowHierarchy: [DeutschlandNodes].[Kreise71Web].[CountyKey71]
       - ColumnHierarchy: [ReportingDate].[Week].[Week] (weekly Meldewoche, 2014--2019)
       - Measure: [Measures].[FallCount_71_Web] (all notified cases; no ReferenceDefinition restriction)
       - Population denominator: Count / Incidence * 1e5 from 2017 [Measures].[FallCount_71_Web]
         and [Measures].[Inzidenz_71_Web], with nearest-year RKI incidence denominator fallback
         (2015/2018) for the three counties reporting zero varicella cases in 2017 (09188, 09763, 16069).
         Note: SurvStat@RKI 2.0 public web service exposes only ReportingDate (not symptom-onset date).
  2. germany_school_holidays_2014_2019.json
     Source: OpenHolidays API (https://openholidaysapi.org/SchoolHolidays), DE-BW .. DE-TH, 2014--2019.
  3. counties_ags.geojson
     Source: https://raw.githubusercontent.com/m-ad/geofeatures-ags-germany/master/geojson/counties.json

Usage:
  python experiment/data/fetch_germany.py          # download all three files
  python experiment/data/fetch_germany.py --check  # verify local files match expected checksums
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = Path(__file__).resolve().parent / "raw" / "germany"

SURVSTAT_URL = "https://tools.rki.de/SurvStat/SurvStatWebService.svc"
OPENHOLIDAYS_URL = "https://openholidaysapi.org/SchoolHolidays"
GEOJSON_URL = "https://raw.githubusercontent.com/m-ad/geofeatures-ags-germany/master/geojson/counties.json"

STATE_MAP = [
    ("01", "DE-SH", "Schleswig-Holstein"),
    ("02", "DE-HH", "Hamburg"),
    ("03", "DE-NI", "Niedersachsen"),
    ("04", "DE-HB", "Bremen"),
    ("05", "DE-NW", "Nordrhein-Westfalen"),
    ("06", "DE-HE", "Hessen"),
    ("07", "DE-RP", "Rheinland-Pfalz"),
    ("08", "DE-BW", "Baden-Württemberg"),
    ("09", "DE-BY", "Bayern"),
    ("10", "DE-SL", "Saarland"),
    ("11", "DE-BE", "Berlin"),
    ("12", "DE-BB", "Brandenburg"),
    ("13", "DE-MV", "Mecklenburg-Vorpommern"),
    ("14", "DE-SN", "Sachsen"),
    ("15", "DE-ST", "Sachsen-Anhalt"),
    ("16", "DE-TH", "Thüringen"),
]

ARCHIVE_EXTRACTION_DATE = "2026-11-19"
NS_MDX = {"b": "http://schemas.datacontract.org/2004/07/Rki.SurvStat.WebService.Contracts.Mdx"}


def _http_get(url: str, timeout: int = 30) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "QD-GMM-Replication/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _soap_olap(year: int, measure: str, col_hierarchy: str, retries: int = 4) -> bytes:
    env = f'''<?xml version="1.0" encoding="utf-8"?>
<soap:Envelope xmlns:soap="http://www.w3.org/2003/05/soap-envelope"
               xmlns:wsa="http://www.w3.org/2005/08/addressing"
               xmlns:rki="http://tools.rki.de/SurvStat/"
               xmlns:mdx="http://schemas.datacontract.org/2004/07/Rki.SurvStat.WebService.Contracts.Mdx">
  <soap:Header>
    <wsa:Action>http://tools.rki.de/SurvStat/SurvStatWebService/GetOlapResultData</wsa:Action>
    <wsa:To>{SURVSTAT_URL}</wsa:To>
  </soap:Header>
  <soap:Body>
    <rki:GetOlapResultData>
      <rki:request>
        <mdx:ColumnHierarchy>{col_hierarchy}</mdx:ColumnHierarchy>
        <mdx:Cube>SurvStat</mdx:Cube>
        <mdx:HierarchyFilters>
          <mdx:KeyValueOfFilterCollectionKeyFilterMemberCollectionb2rWaiIW>
            <mdx:Key><mdx:DimensionId>[PathogenOut].[KategorieNz]</mdx:DimensionId><mdx:HierarchyId>[PathogenOut].[KategorieNz].[Krankheit DE]</mdx:HierarchyId></mdx:Key>
            <mdx:Value><mdx:string>[PathogenOut].[KategorieNz].[Krankheit DE].&amp;[Windpocken]</mdx:string></mdx:Value>
          </mdx:KeyValueOfFilterCollectionKeyFilterMemberCollectionb2rWaiIW>
          <mdx:KeyValueOfFilterCollectionKeyFilterMemberCollectionb2rWaiIW>
            <mdx:Key><mdx:DimensionId>[ReportingDate].[WeekYear]</mdx:DimensionId><mdx:HierarchyId>[ReportingDate].[WeekYear].[WeekYear]</mdx:HierarchyId></mdx:Key>
            <mdx:Value><mdx:string>[ReportingDate].[WeekYear].&amp;[{year}]</mdx:string></mdx:Value>
          </mdx:KeyValueOfFilterCollectionKeyFilterMemberCollectionb2rWaiIW>
        </mdx:HierarchyFilters>
        <mdx:IncludeNullColumns>true</mdx:IncludeNullColumns>
        <mdx:IncludeNullRows>true</mdx:IncludeNullRows>
        <mdx:IncludeTotalColumn>false</mdx:IncludeTotalColumn>
        <mdx:IncludeTotalRow>false</mdx:IncludeTotalRow>
        <mdx:Language>English</mdx:Language>
        <mdx:Measure>{measure}</mdx:Measure>
        <mdx:RowHierarchy>[DeutschlandNodes].[Kreise71Web].[CountyKey71]</mdx:RowHierarchy>
      </rki:request>
    </rki:GetOlapResultData>
  </soap:Body>
</soap:Envelope>'''
    for att in range(retries):
        try:
            req = urllib.request.Request(
                SURVSTAT_URL,
                data=env.encode("utf-8"),
                headers={"Content-Type": "application/soap+xml; charset=utf-8"},
            )
            with urllib.request.urlopen(req, timeout=45) as r:
                return r.read()
        except Exception:
            if att == retries - 1:
                raise
            time.sleep(2 * (att + 1))


def _parse_olap(xml_bytes: bytes):
    root = ET.fromstring(xml_bytes)
    cols = []
    for c in root.findall(".//b:QueryResultColumn", NS_MDX):
        cid = c.findtext("b:ColumnName", "", NS_MDX) or c.findtext("b:Id", "", NS_MDX)
        ccap = c.findtext("b:Caption", "", NS_MDX)
        cols.append([cid, ccap])
    rows = {}
    for r in root.findall(".//b:QueryResultRow", NS_MDX):
        rid = r.findtext("b:RowName", "", NS_MDX) or r.findtext("b:Id", "", NS_MDX)
        rcap = r.findtext("b:Caption", "", NS_MDX)
        vals = []
        for s in r.findall("b:Values/b:string", NS_MDX):
            txt = (s.text or "").strip()
            vals.append(float(txt.replace(",", ".")) if txt else 0.0)
        if rid:
            rows[rid] = {"caption": rcap, "values": vals}
    return cols, rows


def fetch_rki(out_path: Path):
    import datetime
    print("Downloading RKI SurvStat@RKI 2.0 weekly varicella (2014--2019)...", flush=True)
    years_data = {}
    grand_total = 0
    for yr in range(2014, 2020):
        cols, rows = _parse_olap(_soap_olap(yr, "[Measures].[FallCount_71_Web]", "[ReportingDate].[Week].[Week]"))
        yr_rows = {k: {"caption": v["caption"], "counts": [int(round(x)) for x in v["values"]]} for k, v in rows.items()}
        yr_tot = sum(sum(v["counts"]) for v in yr_rows.values())
        if len(yr_rows) < 400 or yr_tot == 0:
            raise RuntimeError(f"SurvStat OLAP parse returned empty/incomplete year {yr}: n_rows={len(yr_rows)}, total={yr_tot}")
        grand_total += yr_tot
        years_data[str(yr)] = {"cols": cols, "rows": yr_rows}
        print(f"  {yr}: {len(cols)} weeks, {len(yr_rows)} counties, total={yr_tot}", flush=True)

    pop_2017 = {}
    _, c17 = _parse_olap(_soap_olap(2017, "[Measures].[FallCount_71_Web]", "[ReportingDate].[WeekYear].[WeekYear]"))
    _, i17 = _parse_olap(_soap_olap(2017, "[Measures].[Inzidenz_71_Web]", "[ReportingDate].[WeekYear].[WeekYear]"))
    for rid, crec in c17.items():
        cnt = float(crec["values"][0]) if crec["values"] else 0.0
        inc = float(i17[rid]["values"][0]) if (rid in i17 and i17[rid]["values"]) else 0.0
        pop = (cnt / inc * 1e5) if (cnt > 0 and inc > 0) else 0.0
        pop_2017[rid] = {"caption": crec["caption"], "count_2017": cnt, "inc_2017": inc, "pop": pop}

    # Nearest-year official RKI incidence-denominator fallback for zero-count counties in 2017
    missing_rids = [rid for rid, rec in pop_2017.items() if rec["pop"] <= 0 and "Unknown" not in rec["caption"]]
    pop_fallback = {}
    for yr_fb in (2018, 2016, 2019, 2015, 2014):
        if not missing_rids:
            break
        _, cy = _parse_olap(_soap_olap(yr_fb, "[Measures].[FallCount_71_Web]", "[ReportingDate].[WeekYear].[WeekYear]"))
        _, iy = _parse_olap(_soap_olap(yr_fb, "[Measures].[Inzidenz_71_Web]", "[ReportingDate].[WeekYear].[WeekYear]"))
        still = []
        for rid in missing_rids:
            cnt = float(cy[rid]["values"][0]) if (rid in cy and cy[rid]["values"]) else 0.0
            inc = float(iy[rid]["values"][0]) if (rid in iy and iy[rid]["values"]) else 0.0
            if cnt > 0 and inc > 0:
                pop_fallback[rid] = {"caption": pop_2017[rid]["caption"], "source_year": yr_fb, "count": cnt, "inc": inc, "pop": cnt / inc * 1e5}
            else:
                still.append(rid)
        missing_rids = still

    import datetime
    payload = {
        "disease": "[PathogenOut].[KategorieNz].[Krankheit DE].&[Windpocken]",
        "extracted_date": datetime.date.today().isoformat(),
        "archive_baseline_date": ARCHIVE_EXTRACTION_DATE,
        "years": years_data,
        "population_2017": pop_2017,
        "population_fallback": pop_fallback,
    }
    out_path.write_bytes(json.dumps(payload, ensure_ascii=False, indent=2).encode("utf-8"))
    print(f"Wrote {out_path} ({out_path.stat().st_size} bytes, total={grand_total})", flush=True)


def fetch_holidays(out_path: Path):
    print("Downloading OpenHolidays KMK German state school holidays (2014--2019)...", flush=True)
    out = {}
    for code, iso, name in STATE_MAP:
        items = []
        for yr in range(2014, 2020):
            url = f"{OPENHOLIDAYS_URL}?countryIsoCode=DE&languageIsoCode=EN&validFrom={yr}-01-01&validTo={yr}-12-31&subdivisionCode={iso}"
            arr = json.loads(_http_get(url).decode("utf-8"))
            items.extend(arr)
        seen = set()
        uniq = []
        for it in items:
            if it["id"] not in seen:
                seen.add(it["id"])
                uniq.append(it)
        out[code] = {"iso": iso, "name": name, "holidays": uniq}
    out_path.write_bytes(json.dumps(out, ensure_ascii=False, indent=1).encode("utf-8"))
    print(f"Wrote {out_path} ({out_path.stat().st_size} bytes)", flush=True)


def fetch_geojson(out_path: Path):
    print("Downloading 5-digit AGS German county polygons...", flush=True)
    raw = _http_get(GEOJSON_URL, timeout=60).replace(b"\r\n", b"\n")
    out_path.write_bytes(raw)
    print(f"Wrote {out_path} ({out_path.stat().st_size} bytes)", flush=True)


def check_local(strict: bool = True):
    import sys
    sys_path_code = str(ROOT / "experiment" / "code")
    if sys_path_code not in sys.path:
        sys.path.insert(0, sys_path_code)
    import data_prep
    verified = data_prep.verify_checksums()
    for k, v in verified.items():
        print(f"  OK  {k:40s}  {v['bytes']:>9d} bytes  {v['sha256']}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="Verify local raw files against checksums")
    ap.add_argument("--verify-download", action="store_true", help="Fail if freshly downloaded files differ from archived checksums")
    args = ap.parse_args()
    if args.check:
        check_local()
        return
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    fetch_rki(RAW_DIR / "rki_varicella_2014_2019.json")
    fetch_holidays(RAW_DIR / "germany_school_holidays_2014_2019.json")
    fetch_geojson(RAW_DIR / "counties_ags.geojson")
    try:
        check_local()
        print(f"All downloaded files match archived SHA-256 (baseline {ARCHIVE_EXTRACTION_DATE}).")
    except RuntimeError as exc:
        if args.verify_download:
            raise
        print(
            f"NOTE: Fresh upstream fetch differs from {ARCHIVE_EXTRACTION_DATE} archive checksum "
            f"(expected when SurvStat updates retroactive case corrections):\n  {exc}"
        )


if __name__ == "__main__":
    main()
