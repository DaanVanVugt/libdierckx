libdierckx_a_OBJECTS = \
			bispev.$(OBJEXT) \
			clocur.$(OBJEXT) \
			cocosp.$(OBJEXT) \
			concon.$(OBJEXT) \
			concur.$(OBJEXT) \
			cualde.$(OBJEXT) \
			curev.$(OBJEXT) \
			curfit.$(OBJEXT) \
			dblint.$(OBJEXT) \
			evapol.$(OBJEXT) \
			fourco.$(OBJEXT) \
			fpader.$(OBJEXT) \
			fpadno.$(OBJEXT) \
			fpadpo.$(OBJEXT) \
			fpback.$(OBJEXT) \
			fpbacp.$(OBJEXT) \
			fpbfout.$(OBJEXT) \
			fpbisp.$(OBJEXT) \
			fpbspl.$(OBJEXT) \
			fpchec.$(OBJEXT) \
			fpched.$(OBJEXT) \
			fpchep.$(OBJEXT) \
			fpclos.$(OBJEXT) \
			fpcoco.$(OBJEXT) \
			fpcons.$(OBJEXT) \
			fpcosp.$(OBJEXT) \
			fpcsin.$(OBJEXT) \
			fpcurf.$(OBJEXT) \
			fpcuro.$(OBJEXT) \
			fpcyt1.$(OBJEXT) \
			fpcyt2.$(OBJEXT) \
			fpdeno.$(OBJEXT) \
			fpdisc.$(OBJEXT) \
			fpfrno.$(OBJEXT) \
			fpgivs.$(OBJEXT) \
			fpgrdi.$(OBJEXT) \
			fpgrpa.$(OBJEXT) \
			fpgrre.$(OBJEXT) \
			fpgrsp.$(OBJEXT) \
			fpinst.$(OBJEXT) \
			fpintb.$(OBJEXT) \
			fpknot.$(OBJEXT) \
			fpopdi.$(OBJEXT) \
			fpopsp.$(OBJEXT) \
			fporde.$(OBJEXT) \
			fppara.$(OBJEXT) \
			fppasu.$(OBJEXT) \
			fpperi.$(OBJEXT) \
			fppocu.$(OBJEXT) \
			fppogr.$(OBJEXT) \
			fppola.$(OBJEXT) \
			fprank.$(OBJEXT) \
			fprati.$(OBJEXT) \
			fpregr.$(OBJEXT) \
			fprota.$(OBJEXT) \
			fprppo.$(OBJEXT) \
			fprpsp.$(OBJEXT) \
			fpseno.$(OBJEXT) \
			fpspgr.$(OBJEXT) \
			fpsphe.$(OBJEXT) \
			fpsuev.$(OBJEXT) \
			fpsurf.$(OBJEXT) \
			fpsysy.$(OBJEXT) \
			fptrnp.$(OBJEXT) \
			fptrpe.$(OBJEXT) \
			insert.$(OBJEXT) \
			parcur.$(OBJEXT) \
			parder.$(OBJEXT) \
			parsur.$(OBJEXT) \
			percur.$(OBJEXT) \
			pogrid.$(OBJEXT) \
			polar.$(OBJEXT) \
			profil.$(OBJEXT) \
			regrid.$(OBJEXT) \
			spalde.$(OBJEXT) \
			spgrid.$(OBJEXT) \
			sphere.$(OBJEXT) \
			splder.$(OBJEXT) \
			splev.$(OBJEXT) \
			splint.$(OBJEXT) \
			sproot.$(OBJEXT) \
			surev.$(OBJEXT) \
			surfit.$(OBJEXT)
libdierckx_a_SOURCES = bispev.f clocur.f cocosp.f concon.f  concur.f cualde.f  curev.f curfit.f dblint.f evapol.f  fourco.f fpader.f  fpadno.f fpadpo.f fpback.f fpbacp.f fpbfout.f fpbisp.f  fpbspl.f fpchec.f fpched.f fpchep.f  fpclos.f fpcoco.f  fpcons.f fpcosp.f fpcsin.f fpcurf.f  fpcuro.f fpcyt1.f  fpcyt2.f fpdeno.f fpdisc.f fpfrno.f  fpgivs.f fpgrdi.f  fpgrpa.f fpgrre.f fpgrsp.f fpinst.f  fpintb.f fpknot.f  fpopdi.f fpopsp.f fporde.f fppara.f  fppasu.f fpperi.f  fppocu.f fppogr.f fppola.f fprank.f  fprati.f fpregr.f  fprota.f fprppo.f fprpsp.f fpseno.f  fpspgr.f fpsphe.f  fpsuev.f fpsurf.f fpsysy.f fptrnp.f  fptrpe.f insert.f  parcur.f parder.f parsur.f percur.f  pogrid.f  polar.f  profil.f regrid.f spalde.f spgrid.f  sphere.f splder.f  splev.f splint.f sproot.f  surev.f  surfit.f
CC = g++
FC = gfortran
OBJEXT = o
RANLIB = ranlib
FCFLAGS = -fdefault-real-8 -fdefault-double-8 -g -Wall

all: libdierckx.a 

libdierckx.a: Makefile $(libdierckx_a_OBJECTS) $(libdierckx_a_SOURCES)
	-rm -f libdierckx.a
	$(AR) cru libdierckx.a $(libdierckx_a_OBJECTS)
	$(RANLIB) libdierckx.a

%.o: %.f Makefile
	$(FC) $< $(FCFLAGS) -c -o $@

clean:
	rm -rf *.o

