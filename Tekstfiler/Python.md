# Torsdag

[NumPy Cheatsheet](/Users/jqc305/Library/Mobile Documents/3L68KQB4HG~com~readdle~CommonDocuments/Documents/Books/Python/Python cheatsheet.pdf)

Gem ting ind i memory. Gem ind i variabel.

I pyhton altid:

for index in files ||| file in files
Altid ental i flertal



lister fylder ikke på samme måde som arrays, arrays er stive og fylder

### Vigtigt

Vi bruger c-python, laves om til byte code, som laves til c, som laves til ML -> virker. Fler mulighed for at fucke op og det tager længere tid pga. flere konverteringer. Interpreters kan konvertere til andre sprog end c (java f.eks). IPython -> Python -> etc. Øger effektivitet.

apply i R analyserer funktion før benytte og optimerer

python kan ikke tåle mix af tab og " ", de fleste editors konvertererer automatisk.

### Functions

Code blocks defineres udelukkende med indents.







Jupyter lab og notebook. Notebook er den gamle, Jupyterlab istedet.

S - Matlab lign. 

Turtles kan kun bevæge sig rundt



## Fredag

Ju()py(python)er® Kan sende data frem og tilbage til hinanden. Bruger mange af de samme libraries (f.eks. forchan). 

del kan fjerne kolonner

brug conda Install

brug argparse pakke til at importere argumenter ind i fil.

```Input til script funktion
import argparse 
parser = argParser(description='tekst')

parser.add_argument('navn', help='beskrivelse af hvad argumentet er')

args= parter.parse_args()

print('whatever du vil sige '+ args.name)  #struktur som med argumenter i nodejs
```

PyPi har alle pakker til download m. pip



[] kan have flere funktioner x [] er lav liste af, x[y], udtag y fra x

Index kan både være 1,2,3 og det svarer til -3,-2,-1

Slice er [ og ) relevant for defence posting

[x:y:z] ->fra x, til y, med afstand z

::2, fra start til slut, men kun hver anden.

[:::-1], baglæns

Arrays er altid firkanter.

List med lister [[],[]] samme som et array, så længe len(inderliste1)=len(inderliste2)



## Notes

* If-else er
  if:
  elif:
  else:
* // floor division, dvs. division hvor den runder af til mindste int og 
* %timeit udregne hvor lang tid noget tager[^1]
* Skriv var? for at få oplysninger om var
* %a=For positive heltal er modulo resten ved [division](https://da.wikipedia.org/wiki/Division_(matematik)). Dvs "overskudet" efter division med int resultat. 27%5=2 fordi 5 går op i 25 og 2 er til overs.
* Python bruger echo så 
  echo x > fil.txt, det er med overwrite. >> er append
  Newline er echo -e og så \n som newline karakteren
* ! kan bruges gentagende gange så "a!=b" er er ikke lig, mens "not a!=b" er lig.
* Indent fungerer istedet for {}
* Start ved 0 og aldrig <=, 0 svarer til =
* Brug func til at undgå at gøre det samme 1 milliard gange
* . operator (så x.append() f.eks.) "calling the function" (ligesom JS), kan bruges uden input. Kan tilføjes til alt (så længe den korrekte inputtype).
  list.sort(my_list) -> 
* List er der intet krav om type, det er der i arrays
* Function? giver help, function?? giver koden
* Hvis man operator på et array, så er det ind i hver felt [1,3,4]*2=[2,6,8]
* Print bruger print.method (giver hvordan det ser ud.)
* Put imports øverst for sig selv!
* ; ved plt.plots fjerner teksten ellers brug plt.show som en voksen person, så det også virker udenfor notebook (kaldes automatisk ved cellestop i jupyter)
* astype('numeric') samme funktion som asNumeric i R
* np.where findes, kan som R, kan bruges til billeder, find alle pladser hvor Lum > x, konverterer det til binært.



[^1]: This module provides a simple way to time small bits of Python code. It has both a Command-Line Interface as well as a callable one. It avoids a number of common traps for measuring execution times. See also Tim Peters’ introduction to the “Algorithms” chapter in the Python Cookbook, published by O’Reilly.