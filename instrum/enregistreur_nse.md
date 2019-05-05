# A trier, source:

si:

```
localuser@cognac:~$ ssh root@192.168.20.204
Unable to negotiate with 192.168.20.204 port 22: no matching key exchange method found. Their offer: diffie-hellman-group1-sha1
```

rajouter dans ssh/.config:
Host 192.168.20.204 
    KexAlgorithms +diffie-hellman-group1-sha1

Pour allumer la source sans émettre de son, il faut ne pas avoir de répertoire `/media/card/sequence/`


# hdrophone

Hydrophone numérique très large bande 1Hz-200kHz

Pas anodin: pas de passe haut.
Si état de mer, les mouvements de la houle seront sur le signal. Pas un pb si le signal n'est saturé pas.
Pb si pilonement de la bouée.
(capteur rbr peut indiquer mvt)

couper les sondeurs de navigation

ylg envoit la doc

connexion: alim + ethernet

adresse manuelle
adresse ordi: 192.168.1.2
sous-masque: 255.255.255.0

adresse hydrophone (browser): 192.168.1.20

10h d'autonomie (vérifié)
128GB de capacité de stockage

Operations:
Régler l'horloge: Set Time, prend l'heure du PC

Recharge:
brancher l'ordi pour monitorer la charge de l'insrument

par défaut chaque .wav fait 1min

Settings: Data collection
Frequence échantillognage: Max Waveform Frequency, 
premier chiffre = fréquence max utilisable, 
chiffre entre parenthèses = frequence échantillognage
pour nous: 25kHz (64kS/s)

Pour lancher la manip:
```
Loggin On
Apply
Apply
```

Arret de l'acquisition:
Setting/Data Collection
```
Logging off
Apply
Apply
```

Téléchargement:
Download

Audacity:
Filtre passe haut 5Hz, 48dB/octave
Gain: saturation pas autorisée, par 10dB

pour couper l'hydrophone: mode standby dans Operations
virer l'alim

synchro?
