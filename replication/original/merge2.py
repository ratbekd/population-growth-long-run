import pandas as pd, numpy as np, warnings; warnings.filterwarnings("ignore")
src=open('merge_tests.py').read()
exec(src.split("# ---------- robustness")[0]); exec("import pycountry\n"+src.split("# ---------- WVS x WPP x World Bank ----------")[1].split('print("\\nIDEAL fertility (D017')[0])
import pycountry_convert as pc
def cont(a3):
    try: return pc.country_alpha2_to_continent_code(pc.country_alpha3_to_country_alpha2(a3))
    except Exception: return {'TWN':'AS','XKX':'EU','TLS':'AS'}.get(a3,'OT')
w['cont']=w.iso3.map({k:cont(k) for k in w.iso3.dropna().unique()})
print(w.drop_duplicates('iso3').cont.value_counts().to_dict())
s=w[(w.S002<=4)&(w.D017<=8)]
print("IDEAL fertility, continent FE")
show("national hab + lnGDPpc, wave + continent FE",s,'D017',['nhab','lgdppc'],['S002','cont'])
show("national + local + lnGDPpc, wave + continent FE",s[s.lhab.notna()],'D017',['nhab','lhab','lgdppc'],['S002','cont'])
r=w[w.age.between(45,70)]
print("REALISED fertility, continent FE")
show("national + local + lnGDPpc, wave + continent FE",r[r.lhab.notna()],'X011',['nhab','lhab','lgdppc'],['S002','cont'])
