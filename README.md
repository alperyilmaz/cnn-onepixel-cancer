# Setup

These steps are needed to run only once.

* `cd raw-data` and then `bash urls`
* `cd misc` and then `docker build . -t alperyilmaz/one-pixel-attack:gpu`

Before each run of `make` you need to source the shell functions defined in `scripts/functions.sh` by running the command `source scripts/functions.sh`

After that, you are ready to go. Just run `make` 