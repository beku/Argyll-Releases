This directory contains the library for sending
rasters to the Google ChromeCast.

Hierarchy:

	cccast.c		Top level actions & receive thread
	ccmes.c			Message handling, uses probuf
	chan/			protobuf encoding
	ccpacket.c		socket write/read

	ccmdns.c		MDNS sign on
	axTLS			
