#
# Regular cron jobs for the xanalysis package
#
0 4	* * *	root	[ -x /usr/bin/xanalysis_maintenance ] && /usr/bin/xanalysis_maintenance
