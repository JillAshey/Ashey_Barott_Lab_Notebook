---
layout: post
title: Fragrammetry models in Metashape
date: '2025-08-26'
categories: Protocols
tags: [Protocols, KBay pairs]
---

## Fragrammetry models in Metashape - surface area extractions 

For the NSF CAREER KBay project, the fragmentss that we use for the PI curves are then photographed by [Spencer Miller](https://www.coralresiliencelab.com/spencer-miller) at CRL using XXXXXXX. We then use these 3D models to extract surface area data from the fragments. 

This post contains information about processing the frag models in Metashape to extract surface area measurements. These instructions were written by Marcelina. 

### Editing the fragrammetry models in Metashape 

The models and images are initially processed by Spencer and then sent to us. Because they are so large, they are stored on the hard drive. 

- To begin processing the models, plug the hard drive into the lab computer and look for the folder of the timepoint of interest. The folder should contain individual folders with images of each frag, as well as the Metashape project file. 
- Open the Metashape project file. This will automatically open Metashape on the computer, as well as a Metashape terminal application. Keep the terminal application open while you work in Metashape. 
- On the left side, there will be the folders with the images of each frag. Click on one of the folders. This should open up the full model of that frag in Metashape. It will look something like this: 

![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/model_pre_edit.png)

- Using the mouse, zoom in/out and spin the model around by clicking and dragging. The sphere is the orientation sphere. You can move the model on an XYZ plane. If you hover over one the lines (red, blue, green) it’ll light up. If you click a specific colored line, you’ll be able to move the model in that specific orientation. 
- On the tool bar, click Model > Free-form selection. The sphere will disappear and a cursor will appear. You can now draw around anything 
	- You can toggle between the selection tool and the orientation sphere by pressing the space bar -- very handy! 
- Draw around the globs (aka model artifacts) that are typically above the frag itself. These should turn pink. Press delete on the keyboard to remove the globs from the model. 
	- If you select something by mistake, click outside the selected part and it will deselect. 

![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/model_glob.png)

- Once the globs are deleted, cut off the base. Don't select all the way up to the frag (see photo below) to avoid cutting off a piece of the coral. 

![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/model_base_selected.png)

- Remove any white pieces left behind. Use the zoom and orientation sphere to move the model around as you select small pieces and delete. This will require a lot of tilting and turning for each frag. 
- Continue to go around the base of the coral until you’re satisfied that you just have the coral

![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/model_editing.png)

- Once the white pieces are removed, turn the model upside down and look inside. Sometimes there are artifacts inside the model that can affect the measurements. 
	- If it is a small glob, leave it 
	- If it is a large glob, remove it. To do this, click Model > Visible selection. This will allow you to select pieces of interest only in the foreground. Change back to free-form selection afterwards. 

![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/frag_hole.png)

- After the inside is clear, fill the model!
	- Press "H" on the keyboard 
	- A popup will appear and the “hole” will be filled in with pink
	- The “level” on the close holes should be 100% 
	- Press “OK” 
	- The popup box will go away but the hole will remain pink. Click outside of the model to solidify it. Once it is solidified, the hole will turn a difference color (blue, gray, etc). 

![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/frag_hole_filled.png)

- Zoom around to double check that everything is properly removed
- Press "Ctrl-S" on the keyboard to save 
- Repeat this process for all frags 

### Measuring the surface area 

- After all frags have been processed as detailed above, the surface area can be measured. 
- Click AddTools > Mesh > Measure 
- This will open up a box. Select the following: 
	- Measure all chunks 
	- Measure = select all 
	- Units = cm
	- Export as CSV
- Click OK 
- The CSV will be exported to the folder where the metashape project lives. 
- Double check the CSV to make sure the values look accurate 



