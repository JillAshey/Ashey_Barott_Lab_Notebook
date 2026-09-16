---
layout: post
title: mtORF amplification 
date: '2026-09-15'
categories: Species ID 
tags: [Heron, Pocillopora]
projects: Heron
---

## mtORF amplification and digestion protocl 

This protocol amplifies the mtORF (mitochondrial open reading frame) of Pocillopora coral DNA using PCR (polymerase chain reaction). The PCR product is then cut with a restriction enzyme to identify Pdamicornis samples. Concept and protocol from [Johnston et al. 2018](https://peerj.com/articles/4355/). Using [Zoe](https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/PCR-Protocol/) and [Hollie](https://github.com/zdellaert/ZD_Putnam_Lab_Notebook/blob/master/protocols/SpeciesID-via-PCR-Sanger-Sequencing.md) protocols as reference. 

#### mtORF Primers

FatP6.1 5′-TTTGGGSATTCGTTTAGCAG-3′
RORF 5′-SCCAATATGTTAAACASCATGTCA-3′

#### mtORF PCR profile
1 cycle: 94 °C for 60 s
30-40 cycles: 94 °C for 30 s, 53 °C for 30 s, and 72 °C for 75 s
1 cycle: 72 °C for 5 min.

### Materials 

- Primers from IDT (see below), stored in -20ºC 
- Ultrapure free water/PCR grade water, stored in -20ºC
- EmeraldAmp® GT PCR Master Mix Cat (# RR310A RR310B) OR BioMix, stored in -20ºC 
- Tsp45I restriction enzyme (NEB), stored in -20ºC 
- 10x rCutSmart buffer (NEB), stored in -20ºC 
- 1kb ladder 
- Loading dye (if needed)
- Thermocyler 
- Strip tubes 
- Pipettes + tips 
- Gloves 

**Important**: Make sure to sterilize all working surfaces with bleach and ethanol before starting to prevent contamination. 

**Important**: All references to water/H20 means PCR-grade, nucelic acid, DNAse, and RNAse free water. 

### Protocol 

#### Primer hydration 

- The primers are shipped dried and so need to be resuspended in liquid before use 
- IDT has a Resuspension Calculation that can be used to hydrate the primers to a concentration of 200 uM. Make a free IDT [account](https://www.idtdna.com/site/account/login?returnurl=%2Fcalc%2Fresuspension) and use the information on the tube to determine how much water to add. 
- Spin down the tube with dehydrated primers before adding water

#### Primer dilution 

- For use in this PCR assay, dilute the hydrated primers from 200 uM to 10 uM. 
- Make a stock solution of 10 uM primer
	- 25 uL of 200 uM primer + 475 uL water 
- Store in -20°C 

#### PCR amplification 

The following formula will be used for the master mix (total volume is 25uL): 

- 10.8 uL water 
- 0.32 uL 10 uM working stock of forward primer 
- 0.32 uL 10 uM working stock of reverse primer 
- 12.55 uL Emerald Master Mix 
- 1 uL of template DNA 

For each amplification, make a positive and negative control. The positive control is a sample of DNA that has been successfully amplified in the past (okay to skip this if you don't have a positive control). The negative control is a sample with no DNA added. The negative control is most important because it will show signs of contamination. 

- Each sample will be amplified once in a 25uL reaction volume 
- Calculate amount of master mix needed. Add number of samples and controls plus some error (5-10%) to get *n* number
	- For example, for 15 samples, make enough for n = 18 
- Make a master mix in a 1.5mL or 5mL tube with these components calculated with your *n* number: 
	- 10.80µL water * *n* =
	- 0.32µL 10µM working stock FatP6.1 primer * *n* =
	- 0.32µL 10µM working stock RORF primer * *n* =
	- 12.55µL Emerald PCR master mix * *n* =
- Vortex and spin down master mix and keep on ice 
- Into labeled PCR tubes, add 24uL of the master mix into each tube
- Add 1uL of sample DNA into each tube. Add 1uL of water for the negative control tube
- Spin down strip tubes 
- Turn on thermocycler and navigate to program for primer set 
	- POC mtORF program: 
		- 1 cycle: 94 °C for 60 s
		- 94 °C for 30 s, 53 °C for 30 s, and 72 °C for 75 s **for 30 cycles **
		- 1 cycle: 72 °C for 5 min.
		- 4°C hold 
- Once program is done, take tubes out and store at 4°C or move to next step 

#### PCR product gel check 

Make a 1% gel to visualize PCR product and successful amplification. 

- Add 75mL new 1x TAE buffer + 0.75g agarose (for small gel) to an Erlenmeyer flask and microwave for ~60 seconds 
	- Keep an eye on the flask in the microwave, as it could boil over during that time. Every ~20 seconds, open the microwave and swirl the flask to mix. Use heat resistant gloves because the flask gets hot. 
- Add 5uL of gel red or 1uL of gel green to the flask
- Set up gel cast mold 
- Pour the gel into the mold and put combs into tray. Let the gel harden/cool until opaque
- Once gel is cool, orient the gel so the comb side is at the top of the box. Take the combs out
- Pour enough "used" TAE buffer into gel box to cover gel with thin layer of liquid 
- Load 4uL of PCR product into the gel (no loading dye should be needed if using the Emerald master mix, as this has loading dye in it already)
- Load 1kb ladder 
- Run gel at 80V for 30 mins
- Visualize gel on gel machine

For POC mtORF, there should be 1 band at ~1000bp. Other bands are signs of potential comtamination and/or off target amplification. 

#### Restriction enzyme digestion 

In the mtORF region, Pdamicornis has a fixed SNP (cytosine; 534 bp), while all other POC had adenine in this position. Therefore, to differentiate Pdamicornis from other POC species, we need to cut the mtORF at this SNP. 

The following formula will be used for the digestion mix (total volume is 10uL): 

- 3.5 uL water 
- 1 uL of 10x rCutSmart buffer 
- 0.5uL of Tsp45I enzyme 
- 5 uL of unpurified PCR product 

- Make a digestion mix in a 1.5mL or 5mL tube with these components calculated with your *n* number: 
	- 3.5µL water * *n* =
	- 1 uL 10x rCutSmart buffer * *n* =
	- 0.5uL of Tsp45I enzyme * *n* =
- Vortext and spin down digestion mix and keep on rice 
- Into labeled PCR tubes, add 5uL of the digestion mix 
- Add 5uL of unpurified PCR product to each tube (the rest of the PCR product can be stored at -20°C)
- Spin down tubes 
- Put in thermocycler and nagivate to digestion program: 65°C for 1 hour, 80°C for ~20 mins, 4°C hold
- Once program is done, run a gel and store the digested product at -20°C

#### PCR gel check 

Run a 2% gel (smaller pore size, good for restriction fragments) for 1 hour at 70V. Load ~4uL from each tube, along with the 1kb ladder. 

If the sample is Pdamicornis, there will be a ~550 bp band (with or without a residual ~1000 bp band) or two bands at 534 bp and ~460 bp. If there is no ~550 bp band, the sample is not Pdamicornis and is some other species. 

### Troubleshooting 

Lots of things can go wrong during PCRs, here are some potential issues that may arise with troubleshooting solutions. 

- Problem: no bands visible on check gel 
	- Issue: low primer concentration
		- Fix: Increase primer volumes from 0.32 to 0.5-0.75uL and adjust water accordingly for a 25uL reaction 
	- Issue: Taq polymerase inactivation 
		- Fix: Aliquot master mix to avoid multiple freeze thaw cycles 
	- Issue: Low DNA input 
		- Fix: Check DNA concentration on Nanodrop. At least 5ng of DNA is needed
- Problem: negative control has bands 
	- Issue: contaminate DNA was introduced to the water, master mix and/or primers 
		- Fix: Throw out working dilutions and stocks and re-make everything after cleaning work station
- Problem: no cutting in Pdam samples 
	- Issue: Enzyme inactivation and/or inhibition. Excess master mix components can interfere with digestion
		- Fix: Keep enzyme as cold as possible while making the digestion mix. 
		- Fix: Ensure PCR product doesn't exceed 50% of total digest volume. Can also perform an ethanol clean up on PCR product before digestion. 
		- Fix: perform sanger sequencing to confirm
- Problem: partial digestion (multiple bands)
	- Issue: insufficient digest time / enzyme - high DNA concentrations and/or low enzyme concentrations can lead to partial digestion 
		- Fix: Extend digestion time to 1.5-2 hours or add additional 0.25uL of enzyme 
		- Fix: perform sanger sequencing to confirm
