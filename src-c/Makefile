LIBNAME=timeslab
LD=ld

OBJS=	timeslab.o

timeslab.so:	$(OBJS)
	@$(LD) $(SHLIBLDFLAGS) -o $(LIBNAME).so $(OBJS)

clean:
	@rm -f *.o *.so *~

realclean:
	@rm -f Makefile *.o *.so
