* make && ~/.local/stata13/stata -b do bench && tail bench.log -n5

set rmsg on
cd ~/bulk/lib/stata-regressby
use /tmp/zz, clear
regressby y x, by(g) bench plugin
use /tmp/zz, clear
regressby y x, by(g) bench
use /tmp/zz, clear
by g: asreg y x, se

* cap net uninstall regressby
* net install regressby, from(https://raw.githubusercontent.com/mcaceresb/stata-regressby/master/) replace
* cd ~/bulk/lib/stata-regressby
*
* * Set up
* clear all
* set obs 10000000
* set seed 123
*
* * Generate a dataset
* gen g = ceil(runiform()*1000)
* gen x = runiform()
* gen y = g + g*x + rnormal()
* sort g
* tempfile t1
* save `t1'
*
* * Test with rangestat
* use `t1', clear
* timer on 1
* * rangestat (reg) y x, interval(g 0 0) by(g)
* timer off 1
* list in 1
*
* * Test with regressby
* use `t1', clear
* timer on 2
* regressby y x, by(g) bench
* timer off 2
* list in 1
*
* * Test with asreg
* use `t1', clear
* timer on 3
* by g: asreg y x, se
* timer off 3
* list in 1
*
* timer list

* cd ~/bulk/lib/stata-regressby
* clear
* set obs 10
* gen x = _n
* gen y = mod(_n, 5)
* gen g = mod(_n, 2)
* gen c = mod(_n, 3)
* sort g c
* l
* preserve
*     regressby y x, by(g) vce(robust) plugin
*     l
* restore
* regressby y x, by(g) vce(robust)
* l
