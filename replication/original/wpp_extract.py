import openpyxl, csv, sys
p="/mnt/user-data/uploads/WPP2024_GEN_F01_DEMOGRAPHIC_INDICATORS_FULL.xlsx"
wb=openpyxl.load_workbook(p, read_only=True, data_only=True)
print("sheets:", wb.sheetnames, flush=True)
ws=wb[wb.sheetnames[0]]
hdr=None; out=csv.writer(open("/home/claude/work2/wpp.csv","w",newline=""))
want=None; n=0
for row in ws.iter_rows(values_only=True):
    if hdr is None:
        cells=[str(c) if c is not None else "" for c in row]
        if any("Total Fertility Rate" in c for c in cells):
            hdr=cells
            def idx(k): return next(i for i,c in enumerate(hdr) if k in c)
            want=[idx("Region, subregion"),idx("ISO3"),idx("Type"),idx("Year"),
                  idx("Total Population, as of 1 July"),idx("Population Density"),idx("Total Fertility Rate")]
            print("header found:", [hdr[i][:40] for i in want], flush=True)
            out.writerow(["name","iso3","type","year","pop_k","density","tfr"])
        continue
    r=[row[i] for i in want]
    if r[2]=="Country/Area": out.writerow(r); n+=1
print("rows written:", n, flush=True)
