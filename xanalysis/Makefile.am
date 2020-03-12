
bin_PROGRAMS = xanalysis dump_dbsav testmpi testcch edt_from_surr

edt_from_surr_SOURCES = edt_from_surr.c

dist_bin_SCRIPTS = retrieve_info.tcl xanalysis_old pslinemerge.pl	\
XLS2CSV.pl xanalysis_shift spikes_by_phase.pl surcontrol.sh		\
setup_passwordless_login cch_rates.pl xanalysis_fix_script.pl		\
xanalysis_nohup xanalysis_run_script xanalysis_make_script

EXTRA_DIST = gen_surr_h.pl commonChangeLog make_interface.pl	\
make_params_h.pl interface.txt surrogateChangeLog debian

dist_doc_DATA = README ChangeLog COPYING LICENSE


xanalysis_SOURCES = menus.f calc_and_write_newshift_cccs.f		\
analyze_data.f display_selected_cycles.f RESP_CYCLE_STATS_3_v4.f	\
rdat_universal.f scalcn.f calc_ach_v4.f calcCTH.f			\
calcCCH_newshift_2.f cardsig_delta2.f respsig6_v4.f intvar4_v2a.f	\
plot_spikes_per_cycle.f showCCHs6_newshift.f showCTHs_v3.f		\
acumsum_nomad_2.f picktypetest2.f betai.f betacfkm.f gammln.f		\
new_draw_button.f enter_date.f locate_region.f new_plot.f		\
locate_region1.f locate_region2.f database_write_new.f stats6.f		\
clear_routines.f calc_cardCCH.f read_info_file.f 			\
ps_plot_1.f miscellaneous_subroutines.f serotonergic_discr_func.f	\
mean_and_sd_routines.f read_and_write_pre_files.f			\
print_and_write_routines.f rootsize.c x2000parameter.defs		\
x2002array.defs xsbf.c ps_starbase.c test_draw_button.f			\
compare_and_convert.f read_and_write_DBSAV.f qdttxt_write.f		\
read_and_write_ref_electrode_data.f version_message.c			\
get_info_IDcodes.f mkdtemp.c mkdtemp.h date.c filesize.c		\
xanalysis_surrogate.c xanalysis_surrogate.h gammafit_search.c unpin.c	\
gen_control.c fpu_check.c gammafit_search.h unpin.h gen_control.h	\
fpu_check.h ps_starbase.h gammafit_search_nopin.c unpin_nopin.c		\
util.h x_cth_cch.c params.h interface.h interface.defs util.c		\
copy_tsfs_files.c write_sta_file.c sbparam.h surrogate_control_c.c	\
pair.c minmax.c host_dialog.c doomed_cells.c pair.h host_dialog.h	\
write_persav_file.c

testcch_SOURCES = testcch.c edt.c cch.c pairloop.c surgen.c tools.c	\
edt.h cch.h pairloop.h surgen.h tools.h

testcch_CFLAGS = `mpicc -showme:compile` $(AM_CFLAGS)
testcch_LDFLAGS = `mpicc -showme:link` $(AM_LDFLAGS) -lrt


x_cth_cch.$(OBJEXT): params.h

dump_dbsav_SOURCES = dump_dbsav.f read_and_write_DBSAV.f

testmpi_SOURCES = testmpi.c
testmpi_CFLAGS = `mpicc -showme:compile` $(AM_CFLAGS)
testmpi_LDFLAGS = `mpicc -showme:link` $(AM_LDFLAGS)

#AM_FFLAGS = +e +U77 -I/opt/graphics/starbase/include
#AM_FFLAGS = +e +U77 +O1 -C
#AM_FFLAGS = +e +U77 +O2 -C -g
#AM_FFLAGS = +e +U77 +O2 +Onolimit -g

if F90
    AM_FFLAGS = +U77 +O2 +Onolimit
    AM_CFLAGS = -mgas
    LDADD = -L/opt/graphics/common -L/opt/graphics/common/lib -lXwindow -lhpgfx -lXhp11 -lX11 $(FLIBS) @GCCLIB@ \
	-L/usr/local/lib -l:libJudy.a -l:libgsl.a -l:libgslcblas.a -lm
    AM_LDFLAGS = -Wl,-N
else
    LDADD =-L/opt/graphics/common -L/opt/graphics/common/lib -lX11 -lm $(FLIBS) @GCCLIB@ -lJudy -lgsl -lgslcblas -lmpi
    AM_FFLAGS = -O2 -fno-underscoring -Wall -fcray-pointer -frecord-marker=4 -fconvert=big-endian -Wno-conversion -fcheck=bounds -ffixed-line-length-0
    AM_CFLAGS = -O2 -Wall --std=c99 -DNOPIN
endif

files: $(xanalysis_SOURCES) $(testcch_SOURCES) $(testmpi_SOURCES) Makefile.am
	ls $(xanalysis_SOURCES) $(testcch_SOURCES) $(testmpi_SOURCES) > files

xanalysis_surrogate.h: gen_surr_h.pl x2000parameter.defs
	./gen_surr_h.pl > xanalysis_surrogate.h

params.h: make_params_h.pl x2000parameter.defs
	./make_params_h.pl

RESP_CYCLE_STATS_3_v4.o: clear_routines.o locate_region2.o miscellaneous_subroutines.o new_draw_button.o new_plot.o
acumsum_nomad_2.o: clear_routines.o locate_region.o miscellaneous_subroutines.o new_draw_button.o print_and_write_routines.o
analyze_data.o: clear_routines.o compare_and_convert.o enter_date.o get_info_IDcodes.o miscellaneous_subroutines.o picktypetest2.o plot_spikes_per_cycle.o print_and_write_routines.o read_and_write_DBSAV.o read_and_write_pre_files.o read_and_write_ref_electrode_data.o read_info_file.o showCCHs6_newshift.o showCTHs_v3.o
betai.o: betacfkm.o gammln.o
calc_and_write_newshift_cccs.o: RESP_CYCLE_STATS_3_v4.o calcCCH_newshift_2.o calcCTH.o calc_ach_v4.o calc_cardCCH.o cardsig_delta2.o clear_routines.o intvar4_v2a.o mean_and_sd_routines.o miscellaneous_subroutines.o new_draw_button.o qdttxt_write.o read_and_write_DBSAV.o respsig6_v4.o scalcn.o serotonergic_discr_func.o showCTHs_v3.o
cardsig_delta2.o: betai.o
compare_and_convert.o: miscellaneous_subroutines.o read_and_write_pre_files.o
convert.o: miscellaneous_subroutines.o read_and_write_pre_files.o
database_write_new.o: miscellaneous_subroutines.o read_and_write_DBSAV.o read_and_write_pre_files.o read_and_write_ref_electrode_data.o
menus.o: analyze_data.o calc_and_write_newshift_cccs.o database_write_new.o enter_date.o miscellaneous_subroutines.o qdttxt_write.o rdat_universal.o read_and_write_DBSAV.o
new_plot.o: miscellaneous_subroutines.o
picktypetest2.o: miscellaneous_subroutines.o
plot_spikes_per_cycle.o: miscellaneous_subroutines.o
print_and_write_routines.o: miscellaneous_subroutines.o ps_plot_1.o
ps_plot_1.o: miscellaneous_subroutines.o
qdttxt_write.o: miscellaneous_subroutines.o read_and_write_DBSAV.o
rdat_universal.o: miscellaneous_subroutines.o
read_and_write_ref_electrode_data.o: miscellaneous_subroutines.o
read_info_file.o: miscellaneous_subroutines.o read_and_write_pre_files.o
respsig6_v4.o: betai.o
serotonergic_discr_func.o: mean_and_sd_routines.o
showCCHs6_newshift.o: acumsum_nomad_2.o clear_routines.o locate_region.o locate_region1.o locate_region2.o miscellaneous_subroutines.o new_draw_button.o new_plot.o plot_spikes_per_cycle.o print_and_write_routines.o stats6.o test_draw_button.o
showCTHs_v3.o: clear_routines.o miscellaneous_subroutines.o new_plot.o
stats6.o: clear_routines.o locate_region.o miscellaneous_subroutines.o new_draw_button.o print_and_write_routines.o
write_qdttxt_intro.o: miscellaneous_subroutines.o


checkin_files = $(edt_from_surr_SOURCES) $(dist_bin_SCRIPTS) $(dist_DOC_DATA) $(EXTRA_DIST) \
					 $(xanalysis_SOURCES) $(testcch_SOURCES) $(dump_dbsav_SOURCES) \
					 $(testmpi_SOURCES) Makefile.am configure.ac ChangeLog


add_head:
	for n in $(checkin_files) ; do \
		touch -r $$n oldtime ; \
		vi $$n ; \
		touch -r oldtime $$n ;\
	done

deb:
	@echo 'Making debian packages'
	make distdir &&\
	cd $(distdir) &&\
	dpkg-buildpackage -Zgzip -uc -us && \
	cd .. && \
	rm -rf $(distdir)

