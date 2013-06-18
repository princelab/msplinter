# msplinter

Predicts how molecules will fragment in a mass spectrometer.  Currently
focused on lipid fragmentation under CID, HCD or PQD.

## Installation

msplinter is built on rubabel, which is built on openbabel.  To get the
necessary openbabel dependencies on ubuntu/debian:

    sudo apt-get install openbabel libopenbabel-dev cmake make curl

Then

    gem install msplinter

See [the "Installing" section](https://github.com/princelab/rubabel) for
complete instructions for installing rubabel and openbabel.

## Commandline Usage

Type 'msplinter' with no arguments to see a help message.

## Copyright

MIT license.  See LICENSE.txt for further details.
