OBJS = soapC.o stdsoap2.o soapClient.o SampleTxbrService.o 
CC = g++
CFLAGS = -Wall -c $(DEBUG)

txbrClient : $(OBJS)
	$(CC) $(OBJS) -o txbrClient 

soapC.o : soapH.h soapC.cpp 
	$(CC) $(CFLAGS)soapC.cpp

stdsoap2.o : stdsoap2.h stdsoap2.cpp
	$(CC) $(CFLAGS) stdsoap2.cpp

soapClient.o : soapH.h soapClient.cpp 
	$(CC) $(CFLAGS) soapClient.cpp 

SampleTxbrService.o : soapH.h SampleTxbrService.cpp txbr.nsmap
	$(CC) $(CFLAGS) SampleTxbrService.cpp

soapH.h : soapStub.h
	cp ../soapH.h .

soapStub.h :
	cp ../soapStub.h .

soapC.cpp :
	cp ../soapC.cpp .

soapClient.cpp :
	cp ../soapClient.cpp .

txbr.nsmap :
	cp ../txbr.nsmap .
	
clean:
	rm *.o txbrClient soapH.h soapC.cpp soapClient.cpp soapStub.h txbr.nsmap

