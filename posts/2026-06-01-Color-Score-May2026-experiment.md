---
layout: post
title: Color score analysis via Fiji/ImageJ
date: '2026-06-01'
projects: NSF CAREER  
---

## Color score analysis of samples from KBay heat stress experiment in May 2026

### Protocol 

#### Naming images 

Before starting on the image analysis, the images of corals need to be renamed so we can keep track of them. The images are stored as jpegs on the Barott Lab Box in this location: `Data/NSF CAREER KBay pairs 2023-2028/Aim 3 - May 2026/Color Score`. Inside the Color Score folder are 7 folders from different days with images in them. The folders are titled YYYYMMDD. 

- Open one of the folders and click on an image. You will see the following: 

xxxx coral image 

- Each picture will have the color card, ruler, ID tag, and the coral on a ceramic plug. 
- Rename each image file so that the name of the image is the ID tag, followed by an underscore and the date. 
	- For example, the image with the coral above would be renamed `PC-242-A1_20260510`. 
	- I already renamed all image files in the `20260506` folder, so use these as a reference for image naming if needed. 
- DO NOT proceed to image analysis before all the images are renamed. 

#### Image analysis 

- Open the color score [spreadsheet](https://docs.google.com/spreadsheets/d/1re59ZWoEaemvJYSRgb9fi3wC84EDLIYtXQd5_irpBmA/edit?gid=0#gid=0) to record data as you go. 
- Download Fiji/ImageJ from [here](https://imagej.net/software/fiji/) onto your local computer. 
- Unzip the folder, and this will create a folder called Fiji. Go into the Fiji folder and there will be an application called Fiji. Click this to open Fiji. 
- You can open images by either going to File > Open and then navigating to the image or by dragging and dropping the image onto the Fiji bar.
- Once the image is up in Fiji, open the histogram analysis option via Analyze > Histogram. 
- Click Live so that it is highlighted in red. This means that you can select different regions of the photo and the histogram will change intensity levels automatically. 
- On the Fiji bar, select the Rectangle tool. 
- Draw a rectangle in as much of the red square as possible. 
- In the histogram window, click through the RGB options and stop at "Red". Record the Mean value in the spreadsheet under "Red.Standard". 
- Repeat the following steps for both the green and blue tape squares. Record the mean value of the "Green" and "Blue" options in their respective columns on the spreadsheet. 
- Select the freehand drawing tool in the Fiji bar (the one that books like a bean). Outline LIVE coral tissue only. Do not select any portions that represent the coral skeleton, the plug, or any part that has a shadow. 
- On the histogram window, click through the RGB options to record the mean red, blue, green AND intensity (unweighted) values. 
- Fill in the rest of the information for the sample in the spreadsheet. 
	- Sample ID - the ID on the tag in each image 
	- Image_Date - the date that the image was taken (corresponds to the date on the folder you are analyzing)
	- Analysis_Date - the date that the image was analyzed
	- Image_Analyzer - name of the person who analyzed the image 
	- Image_name - name of the image file on Box
	- Notes - If there is anything unusual about the image, record it in the notes column

xxx add image of spreadsheet

- Close the image in Fiji and move onto the next one! 
