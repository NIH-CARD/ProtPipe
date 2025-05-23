# README

Install `Apptainer`
```bash
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
```

Build container

```bash
sudo apptainer build protpipe.sif protpipe.def
```