# fonctionnement DRONE ENSTA

## connection au drone

Configuration wifi réseau de l'ordinateur:

```
IPv4: Manual
Adresse IP: 192.168.0.42
Sous-réseau: 255.255.255.0
```

Après avoir brancher le router, tester la bonne connection au flotteur:

```
ping 192.168.0.10
ssh pi@192.168.0.10
[toor]
```

Note: la page de config du router 192.168.0.50 (login=admin, pas de mot de passe), permet de vérfier quels périphériques sont connectés:
`http://192.168.0.50/st_wireless.htm`


## pré-deploiement

- Graisser les parois du piston
- Test piston sinus (si moteur fonction bien): `roslaunch seabot calibration_sinus.launch`
- vérifier état batterie:
`rostopic echo driver/power/battery`
ne doit pas être à 9V, autour de 10V il faut recharger (rapport émilie). 
La batterie est chargée à 12.6V.
- Note: l'environment ros permet l'utilisation de l'autocompletion tab
- Arret propre:  `[ctrl-C]`
(il faut se reconnecter au terminal/screen)
- Arrêt sale: on coupe le courant.
- Si pas de mission pendant 1h le flotteur s'éteint et redémarre et va tuer les terminaux, ne pas être surpris. Pour éviter ça, lancer: `./script_watchdog_alim.sh`
- `seabot/mission/mission_depth_only.xml`: description de la mission
ajuster la date du drone et de la mission: `date -s` est appelé par mission.launch


### recharge de battery

Batterie ont une valeur nominale 11.1V, 5000mAh
Quand on la recharge, il faut indiquer sur le chargeur 3S (3 éléments de 3.1V) et 5A
Batterie Lipo

Attention: brancher pinouille de 4 fils (connecteur de balance) puis +-

Jamais laisser de fil à nue


### Faire le vide ou vérifier que le vide est fait

Pour avoir deux terminaux (avec une session ssh, nécessaire dés fois):

```
screen -a
roslaunch seabot driver.launch
```

`[ctrl-A] [ctrl-D]` détache le terminal

Pour se rattacher:  `screen -r (tab pour autocomplétion)`
(zsh, roslaunch flèche du haut pour dernières commandes)

Liste des topics: `rostopic list`

Liste des noeuds (pour info): `rosnode list`
(attention (ns=)namespace=noeud pas toujours, si préciser ou non dans rospy.publicher(?) )

on vérifie que le capteur interne est présent dans la liste des topic et on lit les données du `rostopic echo /driver/sensor_internal`

pression en hPa, valeur de 1000 au départ, on branche la pompe et on crée une dépression, on descend à 700
(par contre capteur externe est en bar)

pour faire le vide:
dévisser et mettre l'embout de la pompe jusqu'a ce que le joint soit recouvert. 
Faire le vide et monitorer la baisse de pression qui doit être quasi immédiate, descendre légèrement en dessou de 700 (675 par exemple), arreter la pompe et revisser avec l'embout quelques tours avant de finir avec un tournevise classique, jusqu'à atteindre une légère resistance

### définition de la mission

à faire ...

### Test équilibrage

- avoir allumé le drone
- le mettre à l'eau
- faire un `screen -a`
- `roslaunch seabot calibration_equilibrium.launch`
- piston rentre jusqu'à ce que la profondeur dépasse 20cm et s'arrête
- un paramètre à régler: le paramètre de départ du piston (`start_piston_position` dans `seabot/script/calibration_equilibrium.py`), 1400 par exemple.
Avec houle il va falloir peut-être diminuer la valeur
Dans l'idéal il faut faire ça au port avant (i.e. sans houle).

Valeur obtenue à l'issue du test: `Piston Position`.

Cette valeur est à noter et reporter dans `seabot/launch/regulation.launch`, variable `offset_piston`. Mettre un peu moins pour qu'il soit flottant (50)

### Asservisements

Autres paramètres importants pour la régulation:
`K_factor`:  module l'amplitude de la commande,
`K_velocity`: constante de temps (600s approx),
`K_acc`: 0.

(un modèle existe pour trouver les bons paramètres)


### Ouverture du drone

cf photo pour laisser rentrer l'air


## post-déploiement


### récupération

vérifier humidité sur les parois du drone et sur les joints

sudo shutdown now pour éteindre la raspberry manuellement
plus passer l'aimant 2s en face de J1 (photo)

### logs

Les logs se trouvent dans le home: `ls .ros/`

Ceux sont les fichiers .bag (tous les messages) et log des consoles (.ros/log/latest/)

Rappatriement de ces données via scp:

```
scp pi@192.168.0.10:.ros/
```

Pour afficher les données:
`python tools/rosbag/log.py absolutepath/??.bag`
(faire un alias éventuellement)

clique-droit pour zoomer et pan pour explorer les données.


## divers:

consigne=position du piston,
0 (complétement sorti),
2400 (complétement rentré),
à peu prêt (interrupteurs magnétiques ...),
butée de sortie, remise à 0 possible

vérifier date

fichier de lancement principal (en appelle d'autres)
`src/seabot/launch/mission.launch`

- driver.launch: fait fonctionner les capteurs
- filter.launch: processing des données capteurs
- lambert.launch: données gps




## (vieux )à préparer:

- ENSTA

email à lops-tois@listes.ifremer.fr au sujet des recommendations
bouée+ligne pour le drone ENSTA

semaine du 20:
montage + test du proto dédié ifremer
y compris plomb

fin semaine 27, amener drone+accessoires: 
router, chargeurs, pochette, batteries

penser à des missions intéressantes à faire et les "préprogrammer"

script python pour générer depth_only.xml à partir de vecteurs de temps/profondeur éventuellement date


- Ifremer

Emilie montrer recharge batterie à Michel

ensemble bouée+ligne pour drone

graisse pour le piston (quelle graisse - Olivier P.)

préparer une valise pour le flotteur (demander dimensions à Thomas)

tournevis plus boite d'embouts hexagonaux

pompe à vide et tuyau de pompe vide ou deux (olivier)

- indéterminer: plombs d'ajustement 