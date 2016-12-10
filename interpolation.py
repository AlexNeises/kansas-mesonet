from math import pow
from math import sqrt
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import csv
import urllib2
from datetime import datetime, timedelta
import re

def pointValue(x, y, power, smoothing, xv, yv, values):
	nominator = 0
	denominator = 0
	for i in range(0, len(values)):
		dist = sqrt((x - xv[i]) * (x - xv[i]) + (y - yv[i]) * (y - yv[i]) + smoothing * smoothing)
		if dist < 0.0000000001:
			return values[i]
		nominator = nominator + (values[i] / pow(dist, power))
		denominator = denominator + (1 / pow(dist, power))

	if denominator > 0:
		value = nominator / denominator
	else:
		value = -9999
	return value

def invDist(xv, yv, values, xsize = 100, ysize = 100, power = 2, smoothing = 0):
	valuesGrid = np.zeros((ysize, xsize))
	for x in range(0, xsize):
		for y in range(0, ysize):
			valuesGrid[y][x] = pointValue(x, y, power, smoothing, xv, yv, values)
	return valuesGrid

if __name__ == '__main__':
	stations_url = 'http://mesonet.k-state.edu/rest/stationnames/'
	stations_response = urllib2.urlopen(stations_url)
	stations_cr = csv.reader(stations_response)
	station_names = []
	station_latitudes = []
	station_longitudes = []
	for row in stations_cr:
		station_names.append(row[0])
		station_latitudes.append(row[2])
		station_longitudes.append(row[3])
	tm = datetime.now()
	tm = tm - timedelta(minutes=tm.minute % 5, seconds=tm.second, microseconds=tm.microsecond)
	ttest = tm
	stringTime = ttest.strftime('%I:%M %p')
	stringTime = stringTime[1:] if stringTime.startswith('0') else stringTime
	stringTime = stringTime.lower() + ' on ' + ttest.strftime('%B %d, %Y')
	tm = str(tm)
	newtm = re.sub('[^0-9]', '', tm)

	# TEMPORARY
	# newtm = '20161209110500'

	data_url = 'http://mesonet.k-state.edu/rest/stationdata/?stn=all&int=5min&t_start=' + newtm + '&t_end=' + newtm
	data_response = urllib2.urlopen(data_url)
	data_cr = csv.reader(data_response)
	data_station_names = []
	data_weather_data = []
	data_weather_dat2 = []
	keepgoing = True
	for row in data_cr:
		try:
			data_station_names.append(row[1])
			data_weather_data.append(row[11])
			data_weather_dat2.append(row[15])
		except IndexError:
			print 'Data not available, try again in 1 minute.'
			keepgoing = False

	usable_station_data = []

	for x in range(1, len(data_station_names)):
		try:
			new_idx = station_names.index(data_station_names[x])
			usable_station_data.append([station_latitudes[new_idx], station_longitudes[new_idx], data_weather_data[x], data_weather_dat2[x]])
		except IndexError:
			print 'IndexError line 81.'
			keepgoing = False

	if keepgoing:
	
		xv = []
		yv = []
		values = []
		value2 = []

		m1 = interp1d([36, 41], [0, 100])
		m2 = interp1d([-103, -93.6], [0, 100])

		for row in usable_station_data:
			xv.append(m1(float(row[0])))
			yv.append(m2(float(row[1])))
			values.append(float(row[2]))
			value2.append(float(row[3]))

		power = 8
		smoothing = 8

		ti = np.linspace(0, 100, 100)
		XI, YI = np.meshgrid(ti, ti)

		m1fix = interp1d([0, 100], [36, 41])
		m2fix = interp1d([0, 100], [-103, -93.6])
		
		newlat = []
		newlon = []
		newv1 = []
		newv2 = []

		ZI = invDist(xv, yv, values, 100, 100, power, smoothing)
		for (i,j), z in np.ndenumerate(ZI):
			newlat.append(i)
			newlon.append(j)
			newv1.append(z)

		ZI2 = invDist(xv, yv, value2, 100, 100, power, smoothing)
		for (i,j), z in np.ndenumerate(ZI2):
			newv2.append(z)

		with open('winddata.js', 'w') as out:
			out.write('var windData = {timestamp:"%s",x0:-102.0,y0:37,x1:-94.6,y1:40,gridWidth:100.0,gridHeight:100.0,field:[' % (stringTime))
			for idx in range(len(newlat)):
				rad = 4.0*math.atan(1.0)/180.
				u = -newv1[idx]*math.sin(rad*newv2[idx])
				v = -newv1[idx]*math.cos(rad*newv2[idx])
				out.write('%s,%s,' % (u, v))
			out.write(']}')

	# n = plt.Normalize(0.0, 100.0)
	# plt.subplot(111)
	# plt.pcolormesh(YI, XI, ZI, cmap='viridis')
	# plt.colorbar()
	# plt.scatter(yv, xv, 10, values)
	# plt.title('Wind Velocity')
	# plt.xlim(0, 100)
	# plt.ylim(0, 100)

	# plt.show()