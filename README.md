# ARTN
Arizona Robotic Telescope Network

Documents, operating notes, requirements, observing scripts, etc for the Arizona Robotic Telescope Network - 
queue scheduling , remote observing, and automation of UAO telescopes.






[Scott's RTS2 Forck](https://github.com/srswinde/rts2 "Scott's RTS2 fork.")


##Helpfull Info

###Where does RTS2 put data
The path to the image data is defined in the rts2.ini file 
found in /etc/rts2/ directory. Currently you will see 

Currently you will see these definitions in the file:

```
; %b expansion
base_path = "/home/rts2obs/rts2images"

; Images are stored on this path before they are processed.
que_path = "%b/queue/%N/%c/%t/%f"

; The following paths specify where to put images after processing. Every image is first placed
; in queue, and then renamed to target location based on image type, result of astrometry,..
;
; If the following paths are empty/not set, then the images will be not be renamed.
;
; Path for acqusition images
acq_path = "%b/%N/acq/%f"a
; Path for good images - images that were matched with the catalogue
archive_path = "%b/%N/archive/%t/%f"
; Path for bad images - images without match with catalogue.
trash_path = "%b/%N/trash/%t/%f"
; Path for raw skyflats.
flat_path = "%b/%N/skyflats/%f"
; Path for raw darks
dark_path = "%b/%N/darks/%f"
```

For the most part this means that objects in a queue will be stored in 
/home/rts2obs/rts2images/queue/<date>/C0/<object_id>/<unique_name>

