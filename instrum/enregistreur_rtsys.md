
# Matériel requis

- Graisse sylicone pas graisse standard
- cables ethernet et connecteurs aproprié (Ethernet vers USB-C/Thunderbolt)
- alimentation rtsys

Utilisé un code couleur pour les canaux: bleu-blanc-rouge-noir

# Mise en marche

Une pression sur le bouton du boitier de connection mène à 3 cas:

- à une connection à un ordinateur
- à une mission
- ni l'un ni l'autre - s'éteint

Pour se connecter à un ordinateur:

- être connecté en ethernet et déconnecter l'alimentation de préférence

- configuration ethernet ordinateur:

```
IPv4: Manual
Adresse IP: 192.168.20.200
Sous-réseau: 255.255.255.0
```

- dans browser, adresse: `http://192.168.20.204/index` (peut varier d'un enregistreur à un autre)

- si la led est allumé en continue et si la page d'accueil RTSYS apparaît l'ordinateur est connecté. Si la led clignote (1/8-10s), il enregistre.


- dans système configuration:

	- régler l'heure
	- mission schedule: manual, (déactivated en transport)
	- hybrid mode storage
	- active period (temps d'enregistrement): max = 999:59 (3 pression d'affiler pour éteindre)
	- sélectionnées les voix apropriées (faire un petit enregistrement pour tester quelles sont les voix)
	- vérifier la place disponible
	- save & download & shutdown

- 3 pression d'affiler pour éteindre (attendre que le bouton s'éteint à chaque fois)

- 5 heures pour une charge complète

- all in GMT !


# Téléchargement des données

Se connecté à l'instrument et brancher l'alimentation.

Dans filezilla, se connecté à l'adresse:

# Inspection des données:

Via audacity


