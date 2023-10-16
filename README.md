# entanglement-in-geant4
Extension (or rather alteration) of Geant4 to incorporate quantum entanglement, more specifically, Compton scattering of entangled photons.

The provided code is an example of a very simple application, demonstrating the entanglement showing effect during Compton scattering events, i.e., the developed model yielding correct results. It should outline the possible use, and where attention need be paid for any extension to more complex simulations.

Directly addressing the (first) elephant in the room: Originally the code was set up to simulate a Bell's experiment, though a far more appropriate naming concentrates on quantum entangled photons, which was then adopted at some point during this project. Hence half the classes are named something along the lines of 'bell_{usefulNaming}', while the rest are more aptly named. This is a little bit annoying, but was never rectified due to lazyness on the part of yours truly.

Furthermore, the code was created using Geant4 version 10.03.p02, and may not necessarily be compatible other releases. 
