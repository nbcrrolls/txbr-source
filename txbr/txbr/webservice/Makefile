OBJS = soapC.o soapServer.o TxbrWebservice.o stdsoap2.o   
CC = g++
CFLAGS = -Wall -c $(DEBUG)
GENERATEDPATH = ./

txbrWebService : $(OBJS)
	$(CC) $(OBJS) -o txbrWebService 

soapC.o : $(GENERATEDPATH)soapH.h $(GENERATEDPATH)soapC.cpp 
	$(CC) $(CFLAGS) $(GENERATEDPATH)soapC.cpp

stdsoap2.o : $(GENERATEDPATH)stdsoap2.h $(GENERATEDPATH)stdsoap2.cpp
	$(CC) $(CFLAGS) $(GENERATEDPATH)stdsoap2.cpp

soapServer.o : $(GENERATEDPATH)soapStub.h $(GENERATEDPATH)soapServer.cpp 
	$(CC) $(CFLAGS) $(GENERATEDPATH)soapServer.cpp 

soapServerLib.o : $(GENERATEDPATH)soapStub.h $(GENERATEDPATH)soapServerLib.cpp 
	$(CC) $(CFLAGS) $(GENERATEDPATH)soapServerLib.cpp 

TxbrWebservice.o : $(GENERATEDPATH)soapStub.h TxbrWebservice.cpp 
	$(CC) $(CFLAGS) TxbrWebservice.cpp

clean:
	rm *.o *.req.xml *.res.xml *.wsdl *.xsd *.nsmap soapClient* soapC.cpp txbrWebService soapH.h soapObject.h soapServer*.cpp soapProxy.h soapStub.h

