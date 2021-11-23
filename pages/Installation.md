---
layout: page
title: Setting up Jupyter and R kernel in the command line
schemadotorg:
  "@context": http://schema.org/
  "@type": CreativeWork
  "genre": TrainingMaterial
  isPartOf:
      url: "https://gtpb.github.io/IBIP19/"
      name: "IBIP19 - Integrative Biological Interpretation using Proteomics"
---

This should work on a common Ubuntu Linux computer.

Open a terminal with Windows key and then write terminal

Create folder for Training material
`mkdir Training`
`cd Training`

Clone all data onto your computer
`git clone https://github.com/GTPB/IBIP19`

Install jupyter on the computer
`sudo apt install jupyter`

Open R
`R`

Then install the IRkernel package and install it
`install.packages("IRkernel")

Create the links to jupyter
`IRkernel::installspec()`

Get out of R
`q()`
and press "n"


Now the Jupter environment with R kernel should work.

Write
`jupyter notebook`
and press the "new" button to see whether "R" appears


### Back

Back to [main page](../index.md).
