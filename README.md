# hill cipher breaker

finds the hill cipher matrix from a known plaintext and ciphertext pair using qr code alphanumeric charset
then deciphers intercepted messages and ciphers custom messages to send instead

## how it works

* uses 45 symbol set `0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ $%*+-./:`
* detects matrix size automatically (not 1)
* solves linear equations mod 45 to get the hill matrix
* inverts the matrix to decrypt and uses it to encrypt too
