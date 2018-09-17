# Location Privacy Protocol

Implementation of a location protcol based on Spatial Bloom Filters. 
The protocol is designed for a secure multi-party computation setting, where the user and the service provider are mutually distrusting,
 and therefore do not want to disclose private information to the other party. 

## Getting Started

### Prerequisites

This project required you to install the following python-paillier library

```
pip install -r requirements.txt 
```
You will also need the spacial bloom filter class. You can find it here: https://github.com/spatialbloomfilter/libSBF-python

## Bibliography ##
The SBFs protocol have been proposed in the following research papers:
- Luca Calderoni, Paolo Palmieri, Dario Maio: *Location privacy without mutual trust: The spatial Bloom filter.* Computer Communications, vol. 68, pp. 4-16, September 2015. ISSN 0140-3664. [DOI](http://dx.doi.org/10.1016/j.comcom.2015.06.011 "DOI")
- Paolo Palmieri, Luca Calderoni, Dario Maio: *Spatial Bloom Filters: Enabling Privacy in Location-aware Applications*. In: Inscrypt 2014. Lecture Notes in Computer Science, vol. 8957, pp. 16–36, Springer, 2015. [DOI](http://dx.doi.org/10.1007/978-3-319-16745-9_2 "DOI")

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgements 
This project uses several Python modules and in particular: numpy, socket, joblib, collections.

