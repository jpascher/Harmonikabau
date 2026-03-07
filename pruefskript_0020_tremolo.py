#!/usr/bin/env python3
"""
Dok. 0020 — Berechnungsskript
Tremolo: Schwebungsphysik, Typen, Wechselwirkung mit Stimmung

Benötigt: numpy
"""
import numpy as np

f_A4 = 440.0

print("="*80)
print("  Dok. 0020 — Tremolo")
print("="*80)

# 1. Schwebungsfrequenz
print("\n1. SCHWEBUNGSFREQUENZ [Hz] bei Tremolo-Stärke × Grundton")
print("-"*80)
print(f"  {'Δc':>5s} {'Typ':>15s}", end="")
for n,f in [('C3',130.8),('A3',220),('C4',261.6),('A4',440),('C5',523.3),('C6',1047)]:
    print(f" {n:>6s}", end="")
print()
print("-"*80)
for dc,name in [(0.5,'Null'),(1,'sehr leicht'),(2,'leicht'),(3,'mittel'),
                (5,'kräftig'),(8,'Wiener'),(12,'Musette'),(15,'französisch')]:
    print(f"  {dc:4.1f} {name:>15s}", end="")
    for _,f in [('',130.8),('',220),('',261.6),('',440),('',523.3),('',1047)]:
        print(f" {f*dc/1731.23:6.2f}", end="")
    print()

# 2. Wechselwirkung Tremolo × Terzen
print("\n\n2. TREMOLO-SPREIZUNG auf Differenztönen (Terz C4-E4)")
print("-"*60)
f1=261.63; f2=329.63
print(f"  {'Tremolo':>10s} {'Diff min':>10s} {'Diff max':>10s} {'Spreizung':>10s} {'Cent':>8s}")
print("-"*60)
for dc in [0,1,2,3,5,8,12,15]:
    df1=f1*dc/1731.23; df2=f2*dc/1731.23
    diffs=[f2-f1, (f2+df2)-f1, f2-(f1+df1), (f2+df2)-(f1+df1)]
    spread=max(diffs)-min(diffs)
    cent=1200*np.log2(max(diffs)/min(diffs)) if min(diffs)>0 and spread>0 else 0
    print(f"  {dc:8d} Ct {min(diffs):9.2f}  {max(diffs):9.2f}  {spread:9.2f}  {cent:7.1f}")

# 3. Rein vs. temperiert bei Tremolo
print("\n\n3. REIN vs. TEMPERIERT bei verschiedenen Tremolo-Stärken")
print("-"*70)
print(f"  {'Stimmung':>12s} {'Trem.':>6s} {'Schwebung Terz':>16s} {'Gleichmäßig':>12s}")
print("-"*70)
for korr, name in [(0,'Temperiert'),(-7,'Halb-rein'),(-13.7,'Rein')]:
    for dc in [0,2,5,8,12]:
        f2s=f1*2**((400+korr)/1200)
        diff_schweb=abs((f2s-f1)-f1/4)
        df2s=f2s*dc/1731.23; df1s=f1*dc/1731.23
        trem_spr=abs(df2s-df1s)
        gesamt=np.sqrt(diff_schweb**2+trem_spr**2)
        gl="JA" if dc==0 or abs(korr)<1 else "NEIN" if abs(korr)>10 else "teilw."
        print(f"  {name:>12s} {dc:4d} Ct  {gesamt:13.2f} Hz  {gl:>10s}")
    print()

# 4. Tremolo-Frequenz für Begleiter
print("\n4. BEGLEITER: Tremolo-Schwebung bei typischen Akkordtönen")
print("-"*60)
print("   Bass (f₁=130.8 Hz, 3 Cent Tremolo):")
for name,semi in [('C-E gr.Terz',4),('C-G Quinte',7),('C-Eb kl.Terz',3)]:
    f2=130.8*2**(semi/12)
    df1=130.8*3/1731.23; df2=f2*3/1731.23
    d1=f2-130.8; d2=(f2+df2)-(130.8+df1)
    print(f"   {name:>15s}: Diff={d1:.1f} Hz, mit Tremolo={d2:.1f} Hz, Δ={abs(d2-d1):.2f} Hz")

print(f"""

{"="*80}
  ZUSAMMENFASSUNG
{"="*80}

  Schwebungsfrequenz:  f_schweb = f₁ × Δc / 1731,23
  → Proportional zu Grundfrequenz UND Tremolo-Stärke
  → A4 bei 5 Cent: 1,27 Hz, bei 12 Cent: 3,05 Hz

  Tremolo-Typen:
    Null (< 0,5 Cent):  Voraussetzung für Feinstimmung
    Leicht (1–3 Cent):  Steirisch, alpenländisch
    Mittel (3–5 Cent):  Begleiter
    Kräftig (5–8 Cent): Wiener
    Musette (10–15 Cent): Französisch

  Wechselwirkung:
    Tremolo erzeugt eigene Spreizung auf Differenztönen
    Bei > 3 Cent: Reine Stimmung wird verwischt
    Bei > 5 Cent: Ungleichmäßigkeit sofort hörbar
    → Instrumente mit viel Tremolo: TEMPERIERT stimmen
    
  Tremolo steigt proportional mit der Frequenz:
    Gleiche Cent-Verstimmung → höhere Schwebung bei höheren Tönen
    → Tremolo klingt im oberen Diskant schneller als im Bass
    → Manche Stimmer gleichen das aus (weniger Cent oben)
""")


# ══════════════════════════════════════════════════════════════
# 6. TREMOLO-BEZEICHNUNGEN (Industrie)
# ══════════════════════════════════════════════════════════════
print("\n\n6. TREMOLO-BEZEICHNUNGEN (Industrie/Hersteller)")
print("-"*80)
print(f"{'Hz (A4)':>7s} {'Cent':>6s} {'Bezeichnung':>30s} {'Stil/Genre':>30s}")
print("-"*80)
data = [
    (0, 0, 'Unison / Dry / Secco', 'Klassik, Balkan'),
    (0.5, 2, 'Concert / Violin / Swing Secco', 'Jazz, Tango, Cleveland Polka'),
    (1.0, 4, 'Swing / Swing Mosso', 'Gypsy Jazz, Brazilian, Klezmer'),
    (2.0, 7, 'Demi-Swing / Mezzo Swing', 'Irish'),
    (2.5, 10, 'American / Americano', 'Cajun, Québécois'),
    (3.0, 12, 'Moderate / Slovenisch', 'Tex-Mex, Alpin'),
    (4.0, 15, 'Standard / German / Italian', 'Polka, Continental'),
    (5.0, 18, 'Fast / Modern French', 'Modernes Musette'),
    (6.0, 22, 'Very Fast / Old French', 'Altes Musette'),
    (7.0, 25, 'Extremely Fast / Scottish', 'Scottish'),
]
for hz, ct, name, genre in data:
    print(f"  {hz:5.1f}  {ct:5d}  {name:>30s}  {genre:>30s}")

# 7. Akkord-Verstimmung Beweis
print("\n\n7. AKKORD-VERSTIMMUNG durch Tremolo (Cross-Kombinationen)")
print("-"*60)
f1=261.63; f3=f1*5/4; f5=f1*3/2
for dc in [0,3,5,8,12]:
    f3t=f3*2**(dc/1200); f1t=f1*2**(dc/1200)
    terz=1200*np.log2(f3/f1)
    cross1=1200*np.log2(f3t/f1); cross2=1200*np.log2(f3/f1t)
    print(f"  Tremolo {dc:2d} Ct: Terz ±{abs(cross1-terz):.1f} Cent Spreizung = {abs(cross1-cross2):.0f} Cent")

# 8. Tremolo-Verläufe (S-Form Beispiel)
print("\n\n8. S-FORM Tremolo-Verlauf (Beispiel)")
print("-"*50)
# Typisch: unten 1-2 Hz, Mitte 2-3 Hz, oben 4-6 Hz
import numpy as np
noten=['F3','C4','A4','C5','A5','C6']
freqs=[174.6, 261.6, 440.0, 523.3, 880.0, 1047.0]
# S-Form: Cent steigt von 3 unten auf 8 oben
cents_s=[3, 4, 5, 6, 7, 8]
print(f"  {'Ton':>5s} {'f [Hz]':>8s} {'Cent':>6s} {'Hz Schwebung':>14s}")
for n,f,c in zip(noten,freqs,cents_s):
    hz=f*c/1731.23
    print(f"  {n:>5s} {f:7.1f}  {c:5d}   {hz:12.2f}")
