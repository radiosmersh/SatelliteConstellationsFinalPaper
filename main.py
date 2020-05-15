import spacetrack.operators as op
from spacetrack import SpaceTrackClient
from tletools import TLE
from tletools.pandas import load_dataframe
from datetime import datetime
from skyfield.api import Topos, load
from skyfield.units import Distance
from matplotlib import pyplot as plt
from astropy.visualization import hist
from itertools import combinations
from math import degrees, pi


USERNAME = 'REDACTED'
PASSWORD = 'REDACTED'

st = SpaceTrackClient(identity=USERNAME, password=PASSWORD)
data_generator = st.tle_latest(iter_lines=True, ordinal=1, epoch='>now-30', format='3le')

with open('data.txt', 'w') as f:
    for line in data_generator:
        f.write(line + '\n')

df = load_dataframe('data.txt')

df['a'] = df['n'].apply(lambda n: (8681663.653 / n) ** (2/3))
df['ah'] = df['a'] * (1 - df['ecc']) - 6371
df['ph'] = df['a'] * (1 + df['ecc']) - 6371


satellites_tle = []

satellites = load.tle_file('data.txt')

print(len(satellites))

ts = load.timescale()
t = ts.now()

# Высокоэллиптические орбиты
HEO = df[(df['ah'] <= 2000) & (df['ph'] >= 35786)]

# «Молния»
molniya = HEO[(HEO['inc'].between(60.8, 64.8)) & (HEO['n'].between(2.004, 2.01))]

# «Тундра»
tundra = df[(df['argp'].between(240, 300)) & (df['n'].between(1.00173790901, 1.00373790901)) & (df['ecc'].between(0.25, 0.4))]

print('Количество спутников на ВЭО: ', len(HEO.index))
print('Количество спутников c орбитой типа «Молния»: ', len(molniya.index))
print('Количество спутников c орбитой типа «Тундра»: ', len(tundra.index))


# ГСО
GSO = df[(df['inc'].between(0, 2, inclusive=True)) & (df['ecc'].between(0, 0.1, inclusive=True)) & (df['n'].between(0.9, 1.2, inclusive=True))]
inclined_GSO = GSO[GSO['inc'] > 0.15]

print('Количество спутников c орбитой типа «Молния»: ', len(GSO.index))
print('Количество спутников c орбитой типа «Тундра»: ', len(inclined_GSO.index))


by_number = {sat.model.satnum: sat for sat in satellites}

GSO_longitudes = []

GSO_satellites = []

for index, tle_object in GSO.iterrows():
    satellite = by_number[int(tle_object['norad'])]
    GSO_satellites.append(satellite)
    geometry = satellite.at(t)
    subpoint = geometry.subpoint()
    d, m, s = subpoint.longitude.dms()
    GSO_longitudes.append(round(d + m / 60 + s / 3600))

longitudes = {}

for lon in range(-181, 181):
    longitudes[lon] = GSO_longitudes.count(lon)

fig, ax = plt.subplots()
ax.bar(longitudes.keys(), longitudes.values())
ax.set_title(u'Распределение геостационарных КА по долготе')
ax.set_xlabel(u'Долгота,°')
ax.set_ylabel(u'Кол-во спутников')
plt.savefig('Распределение геостационарных КА по долготе.png')

min_distances = []

for s1 in GSO_satellites:
    p1 = s1.at(t)
    min_dist = Distance(au=100).km
    for s2 in GSO_satellites:
        if s2 == s1:
            continue
        p2 = s2.at(t)
        dist = (p2 - p1).distance().km
        if dist < min_dist:
            min_dist = dist
    min_distances.append(min_dist)

fig, ax = plt.subplots()
hist(min_distances, bins='scott', density=True, alpha=0.5, ax=ax)
ax.set_title(u'Распределение КА по расстоянию между «ближайшими соседями»')
ax.set_xlabel(u'Расстояние, км')
ax.set_ylabel(u'Частота')
plt.savefig('Распределение КА по расстоянию.png')


orbit_heights = []

eccos = []

LEO_sat_count = 0
SSO_sat_count = 0

for satellite in satellites:
    height = satellite.at(t).subpoint().elevation.km
    if height <= 2000:
        orbit_heights.append(height)
        eccos.append(satellite.model.ecco)
        LEO_sat_count += 1
        
        if ((height >= 600) and (height <= 800)) and \
            (abs(degrees(satellite.model.inclo) - 98) <= 0.5):
                T = 2 * pi / satellite.model.no_kozai
                if (T >= 96) and (T <= 100):
                    SSO_sat_count += 1

fig, ax = plt.subplots()
hist(orbit_heights, bins='scott', density=False, alpha=0.5, ax=ax)
ax.set_xlim(left=0, right=2000)
ax.set_title(u'Распределение низкоорбитальных КА по высотам орбиты')
ax.set_xlabel(u'Высота орбиты, км')
ax.set_ylabel(u'Частота')
plt.savefig('Распределение низкоорбитальных КА по высотам орбиты.png')

fig, ax = plt.subplots()
ax.set_xlim(left=0, right=0.1)
hist(eccos, bins='scott', density=False, alpha=0.5, ax=ax)
ax.set_title(u'Распределение низкоорбитальных КА по эксцентриситетам')
ax.set_xlabel(u'Высота орбиты, км')
ax.set_ylabel(u'Частота')
plt.savefig('Распределение низкоорбитальных КА по эксцентриситетам.png')

print(LEO_sat_count, SSO_sat_count)