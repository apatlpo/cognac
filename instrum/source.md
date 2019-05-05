# acoustical source, rtsys


## basic usage

Ethernet configuration on computer:
```
Manual IPv4 configuration
IP Address: 192.168.20.200
Subnet Mask: 255.255.255.0
```
Connect to in browser at the following address: [http://192.168.20.204](http://192.168.20.204)

Same with Bullet wifi antenna, but beware:

- do not swap ethernet cables: thicker cable goes to bullet
- make sure wifi has not been turned off in `wifi_config.json` file

Do not charge with the on/off cable plugged.

In order to fire the source without sound emission, you muste not have a repertory `/media/card/sequence/`

## potential issues

If:

```
localuser@cognac:~$ ssh root@192.168.20.204
Unable to negotiate with 192.168.20.204 port 22: no matching key exchange method found. Their offer: diffie-hellman-group1-sha1
```

add in `ssh/.config`:

```
Host 192.168.20.204
    KexAlgorithms +diffie-hellman-group1-sha1
```
