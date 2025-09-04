---
layout: post
title: CoralNet instructions
date: '2025-09-03'
categories: Analysis
tags: [CoralNet, Photos]
projects: NSF CAREER
---

## CoralNet benthic quadrat photo analysis 

[CoralNet](https://coralnet.ucsd.edu/source/) is a resource for benthic image analysis. It allows users to upload images and classify the benthic cover based on labels of their choosing. You can read more about CoralNet in [Beijbom et al. 2015](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0130312) and [Williams et al. 2019](https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2019.00222/full) if interested. 

### Field methods 

Benthic transects are done in Kaneohe Bay, Hawaii at Patch Reef 13 during our seasonal sampling trips. We follow methods described in [Matsuda et al. 2020](https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2020.00178/full). 

- Four surveys are conducted at Reef 13 (two surveys at 3m depth, two surveys at 1m depth). 
- A 40m transect is laid parallel to the reef crest. 
- Every 2m, a photoquadrat is placed on the benthos and a photo is taken of the photoquadrat. 

In total: 4 surveys x 20 pictures per transect = 80 photos per sampling trip 

### Naming photos 

Lots of the photos are already on CoralNet and annotated. The ones that we will work on will be in the Barott lab Box in the following locations: 

- `Data > NSF CAREER KBay pairs 2023-2028 > Benthic photoquad images`
	- `2024-10-08_KBay_Reef13_photoquad_images` - photos are in transect folders, but not named 
	- `2025-01-08_KBay_Reef13_photoquad_images` - photos are not in transect folders and are not named 
	- `2025-05-12_KBay_Reef13_photoquad_images` - photos are in transect folders, but not named 
- `Data > KBay coral pairs and temp time series > Benthic photoquadrat surveys`
	- `2015-10 images` - photos are named (ONLY upload Site 13)
	- `2015-12 images` - photos are named (ONLY upload Site 13)
	- `2016-03 images` - photos are named (ONLY upload Site 13)
	- `2016-09 images` - photos are named (ONLY upload Site 13)
	- `2015-10 images` - photos are named (ONLY upload Site 13)
	- `2017-06 images` - photos are named (ONLY upload Site 13)

Some photos may already be named but if they not, the naming scheme will go like this: `Date_Reef13_XXm_TransectY_ZZ` where `Date` is the date the photo was taken (should be in the folder name), `XX` is the depth at which the photo was taken (1 or 3m), `Y` is the transect number (1, 2, 3, or 4), and `ZZ` is the number of the photo (1-20). 

There may be duplicates of some pictures. Pick the best one out of the duplicates to rename and use in analysis. Additionally, some photos will be in folders based on transect number. If the photos are not in transect folders, you can put "Transect" in the photo name instead of "Transect1" etc. 

### Registering for CoralNet 

- Make an account on CoralNet 
- Send me (Jill) your username and I will add you to the KBayBleach2019 repository

![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/KBay_coralnet.png)

On the KBayBleach2019 main page, an image status box can be seen on the right. This shows you how many images there are total and how many have been annotated (i.e., confirmed). If images have not yet been annotated, they will show up as Unclassified here. 

### Uploading to CoralNet 

After the photos are named, they can be uploaded to CoralNet. 

- Go to the KBayBleach2019 repository and click Upload 
- Click upload images 
- Choose files that you want to upload 
- Click Start upload. This may take a couple of minutes depending on how many files you are uploading. 
- Once the images are uploaded, click Manage image metadata 
- For each image, add the following metadata information: Date (that the image was taken), Site (Reef13), Transect (1-4), Image (1-20), CoralDepth (1 or 3m). 
	- This information should already be in the image name but its good to have the information in multiple places. 

### Annotating the photos 

Once the images are uploaded, it is time to annotate!

- Click on Images 
- Photos that have not yet been analyzed will have a red square around them. Photos with a green square mean that the image has been annotated. 

![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/KBayBleach2019_coralnet_images.png)

- Click on a photo with a red square. This will take you to a page with just the image. 

![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/photo_example.png)

- Scroll down to Annotation and point location status
- Click Edit next to Annotation area
- Adjust the yellow borders to fit completely inside the quadrat
	- Depending on the angle of the photos, some photos may be more difficult and may cut off some area inside the quadrat 

![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/redraw_boundary.png)

- Click Save changes. This will take you back to the image (it won't look any different).  
- Click Annotation Tool at the top right of the page. 
- This will generate 50 random points inside of the newly drawn border.

![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/points.png)

- Click on each point to annotate it (more information on annotations below).
	- Use the zoom in and out buttoms to get closer to the point (upper right hand side of the image).
	- Use the brightness and contrast sliders if needed. 
	- Click on the single point button to see one point at a time. 
- Click through the points to annotate them all. 
- Once all points are annotated, click save! 

Hooray, now that image is annotated and classified! It should now have a green square around it when on the Images page. 

### Annotation and label examples 

There are a lot of potential labels to use and some can be redundant. To avoid confusion and maintain consistency, we will use the labels detailed below. 

| Label code | Name | Functional Group | Example |
|---------|---------|---------|---------|
| B_PorComp | Bleached Porites compressa | Hard coral | ![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/B_PorComp_example.png) |
| C-rubble | Broken coral rubble | Other | ![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/C-rubble_example.png) |
| CCA | CCA (crustose coralline algae) | Algae | ![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/CCA_example.png) |
| D_coral | Dead coral | Hard substrate | ![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/D-coral_example.png) |
| Dictyosph | Dictyosphaeria | Algae | ![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/Dictyosph_example.png) |
| MCAP Pale | Montipora capitata paled | Hard coral | ![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/MCAP_Pale_example.png) |
| MCAP_Pig | Pigmented Montipora capitata | Hard coral | ![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/MCAP_Pig_example.png) |
| Mgrand | Mycale grandis | Other invertebrates | ![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/Mgrand_example.png) |
| Moncap_BL | Montipora capitata bleached | Hard coral | ![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/Moncap_BL_example.png) |
| PCMP Pale | Porites compressa paled | Hard coral | ![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/PCMP_pale_example.png) |
| PCOM_Pig | Pigmented Porites compressa | Hard coral | ![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/PCOM_Pig_example.png) |
| Sand | Sand | Soft substrate | ![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/sand_example.png) |
| SHAD | Shadow | Other| ![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/SHAD_example.png) |
| TT | Transect Tape | Other| ![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/TT_example.png) |
| Turf_algae | Turf algae | Algae| ![](https://raw.githubusercontent.com/JillAshey/Ashey_Barott_Lab_Notebook/refs/heads/main/images/coralnet/Turf_algae_example.png) |

Other labels that may be used could be [PocDam](https://www.marinelifephotography.com/corals/cauliflower/pocillopora-damicornis.htm) (Pocillopora damicornis, now Pocillopora acuta; hard coral), [zoanthid](https://www.waikikiaquarium.org/experience/animal-guide/invertebrates/zoanthids/) (other invertebrates), or [LPurp](https://www.marinelifephotography.com/corals/faviidae/leptastrea-purpurea.htm) (Leptastrea purpurea; hard coral). These will be much less common.

If you can't figure out how to label a point, send it to me or label it as Unk (Unknown). Some of the images or image features may be grainy or difficult to see, but do the best you can!