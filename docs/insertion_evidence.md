# Evidence for cuff apposition and strip insertion

This note records the primary measurements used to decide which insertion scenarios can be drawn. None of these studies measures a penetration probability for the supplied 48-site device in the human median nerve.

## Cuff reshaping

Tyler and Durand tested a flat interface nerve electrode that reshapes peripheral nerves for selective stimulation. This supports representing the cuff as a fixed shape and the nerve as the locally reshaped body. It does not determine the 0.12 mm indentation selected in `quicklook3DParams.m`. That value remains an illustrative input.

Source: Tyler DJ, Durand DM, *IEEE TNSRE* (2002), DOI [10.1109/TNSRE.2002.806840](https://doi.org/10.1109/TNSRE.2002.806840).

## Separation between fascicles

Tyler and Durand's SPINE used blunt elements to separate and rearrange epineurial tissue without compromising the perineurium. The result motivates a geometric path between fascicles. It does not establish whether this device's strip bends, fascicles move, or both.

Source: Tyler DJ, Durand DM, *IEEE TRE* (1997), [PubMed 9086385](https://pubmed.ncbi.nlm.nih.gov/9086385/).

Harreby and colleagues implanted six transverse intrafascicular multichannel electrodes in the median nerves of four Göttingen minipigs. Their histology placed all six between fascicles after 32–37 days. This is a device- and animal-specific observation, not a population estimate of the probability of fascicle penetration for the project device.

Source: Harreby KR et al., *Artificial Organs* (2015), DOI [10.1111/aor.12347](https://doi.org/10.1111/aor.12347).

A 2025 ex vivo study guided a penetrating cuff contact between sciatic nerve fascicles and checked its position with micro-CT. The controlled placement is evidence that an interfascicular route is achievable for that setup. It provides no probability for an unguided strip in a human median nerve.

Source: Ben M'Rad I et al., *Bioelectronic Medicine* (2025), DOI [10.1186/s42234-025-00193-6](https://doi.org/10.1186/s42234-025-00193-6).

## Insertion mechanics and transfer limits

Isolated rabbit sciatic epineurium had a Young's modulus of 0.4 ± 0.1 MPa and perineurium 3.0 ± 0.3 MPa. The study measured insertion forces using sharp and blunt 100 μm steel probes at 2 mm/s. These values refer to isolated rabbit tissue, not to this device, its strip substrate, or an in-situ human nerve. They are therefore stored as literature observations and are not used to set the route in the current code.

Source: *Biomechanical characterization of isolated epineurial and perineurial membranes of rabbit sciatic nerve*, *Journal of Biomechanics* (2022), DOI [10.1016/j.jbiomech.2022.111058](https://doi.org/10.1016/j.jbiomech.2022.111058).

A needle experiment on excised rabbit nerve reported dimpling before insertion and about 2 mN for epineurial puncture with a 50 μm electro-sharpened tungsten needle. Its insertion speed was uncontrolled and its probe differs from this strip. The force does not provide a threshold for project-device fascicle penetration.

Source: Jensen W et al., *Measurement of intrafascicular insertion force of a tungsten needle into peripheral nerve* (2001), [author PDF](https://www.isip.uni-luebeck.de/fileadmin/files/publications/Jensen_ieee2001.pdf).

## Model status

`run_3d_quicklook("deflect")` searches centreline paths around fixed fascicles and rejects paths that overlap their perineurium plus an assumed strip envelope. `run_3d_quicklook("penetrate")` records geometric overlap on a straight centreline. Neither run solves tissue displacement, beam bending, puncture, or insertion probability. A mechanical model requires the strip's flexural rigidity and tip shape, insertion constraints, nerve constitutive parameters, and anatomical geometry at the chosen level. The unspecified values must be measured, supplied by the group, or explicitly explored as sensitivity ranges.

## Direction of deflection and impact angle

A micro-ultrasound study of three human cadaver median nerves observed fascicle displacement and rotation during needle insertion. Four attempts were needed to pass through the perineurium in a fresh sample. This directly supports fascicle movement as one possible response. The instrument was a needle and the specimens were ex vivo, so these observations do not determine the deformation of this project's flexible strip.

Source: *The feasibility of micro-ultrasound as a tool to image peripheral nerves*, *Anaesthesia* (2016), DOI [10.1111/anae.13708](https://doi.org/10.1111/anae.13708).

A pig brachial-plexus micro-ultrasound experiment reported epineurial tenting and nerve rotation. In that preparation, trajectories slightly away from the perpendicular tended to slide off the epineurium. This is evidence about the outer nerve surface, not a measured relation between incidence angle and penetration of an individual fascicle by a flexible strip.

Source: *Real-time visualisation of peripheral nerve trauma during subepineural injection in pig brachial plexus using micro-ultrasound*, *British Journal of Anaesthesia* (2021), [full paper](https://www.bjanaesthesia.org/article/S0007091221002373/pdf).

The current path planner keeps fascicle centres fixed while shifting the strip, because the relative motion cannot be resolved from the available mechanical data. In the real nerve, both the fascicles and strip may move. The planner's chosen bend amplitude is a geometric sensitivity input and must not be interpreted as an observed mechanical displacement.
