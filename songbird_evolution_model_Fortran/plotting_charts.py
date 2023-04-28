

import pandas, numpy
wine_data = pandas.read_csv('Book1.csv', sep = ';', header = None)
wine_data_ = wine_data
wine_data = numpy.array([x.split(';') for x in wine_data_[0]])
wine_data



