from matplotlib import dates
import datetime

data, times, ratios = import_time_data()
new_times = []
for time in times:
    new_time = dates.num2date(dates.epoch2num(float(time) + 1298962800))
    new_times.append(new_time)


hfmt = dates.DateFormatter('%m/%d %H:%M')
fig, ax = plt.subplots()
plt.plot(new_times, ratios, c='r')

fig.subplots_adjust(bottom=0.2)
ax.xaxis.set_major_formatter(hfmt)
fig.autofmt_xdate()


plt.show()


