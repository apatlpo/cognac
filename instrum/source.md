# acoustical source, rtsys


## basic usage

Check at deployment that the GPS data is logged, e.g. `GPS :: $GNRMC` lines in `/media/card/records/mission... .txt`

**Ethernet** configuration on computer:

```
Manual IPv4 configuration
IP Address: 192.168.20.200
Subnet Mask: 255.255.255.0
```

Connect to in browser at the following address: [http://192.168.20.204](http://192.168.20.204)

**wifi**: with bullet antenna, set address with automatic DHCP, but beware:

- do not swap ethernet cables: thicker cable goes to bullet
- make sure wifi has not been turned off in `wifi_config.json` file

Do not charge with the on/off cable plugged.

In order to fire the source without sound emission, you must not have a repertory `/media/card/sequence/`

## potential issues

- the source will shutdown automatically if `/media/card/emission_config.yaml` does not exist and if no ethernet connexion is established.

- Error when connecting to source from a mac via ssh:

```
localuser@cognac:~$ ssh root@192.168.20.204
Unable to negotiate with 192.168.20.204 port 22: no matching key exchange method found. Their offer: diffie-hellman-group1-sha1
```

**Un**comment in `/etc/ssh/ssh_config` the following lines:

```
Ciphers aes128-ctr,aes192-ctr,aes256-ctr,arcfour256,arcfour128,aes128-cbc,3des-cbc MACs hmac-md5,hmac-sha1,umac-64@openssh.com,hmac-ripemd160
```

and **add** the following ones:

```
Host 192.168.20.204
    KexAlgorithms +diffie-hellman-group1-sha1
```
