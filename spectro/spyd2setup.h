#ifndef SPYD2SETUP_H

/* 
 * Argyll Color Correction System
 *
 * ColorVision Spyder 2 related software.
 *
 * Author: Graeme W. Gill
 * Date:   19/10/2006
 *
 * Copyright 2006 - 2007, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* This file is only included in top utilities that need to */
/* be able to access the Spyder 2 colorimeter. This provides */
/* a mechanism for ensuring that only such utilities load the */
/* proprietary Spyder firmware, as well as providing a means to */
/* detect if the spyder driver is going to be funcional. */

extern unsigned int *spyder2_pld_size;			/* in spyd2.c */
extern unsigned char *spyder2_pld_bytes;

/* Return 0 if Spyder 2 firmware is not avaialable */
/* Return 1 if Spyder 2 firmware is available from an external file */
/* Return 2 if Spyder 2 firmware is part of this executable */
int setup_spyd2() {
	unsigned int size, rsize;
	static loaded = 0;		/* We attempted loading from spyd2PLD.bin */
	FILE *fp;
	int i;

	/* Spyder 2 Colorimeter Xilinx XCS05XL firmware pattern. */
	/* This is a placeholder in the distributed files. */
	/* It could be replaced with the actual end users firmware */
	/* by using the spyd2trans utility, but normally the spyd2PLD.bin */
	/* file is loaded instead. */

#include "spyd2PLD.h"

	spyder2_pld_size = &pld_size; 
	spyder2_pld_bytes = pld_bytes; 

	/* If no firmware here, see if there is a binary file to load from. */
	if ((pld_size == 0 || pld_size == 0x11223344) && loaded == 0) {
		char binpath[MAXNAMEL+1];
		loaded = 1;
		
		for (;;) {
			if (strlen(exe_path) + strlen("spyd2PLD.bin") > MAXNAMEL)
				break;				/* oops */
			strcpy(binpath, exe_path);			/* Use global */
			strcat(binpath, "spyd2PLD.bin");
	
			/* open binary file */
#if defined(O_BINARY) || defined(_O_BINARY)
			if ((fp = fopen(binpath,"rb")) == NULL)
#else
			if ((fp = fopen(binpath,"r")) == NULL)
#endif
				break;

			/* Figure out how file it is */
			if (fseek(fp, 0, SEEK_END)) {
				fclose(fp);
				break;
			}
			size = (unsigned long)ftell(fp);

			if (size > pld_space) 
				size = pld_space;
		
			if (fseek(fp, 0, SEEK_SET)) {
				fclose(fp);
				break;
			}
		
			if (fread(pld_bytes, 1, size, fp) != size) {
				fclose(fp);
				break;
			}
			pld_size = size;
//printf("~1 bytes = 0x%x 0x%x 0x%x 0x%x\n",
//pld_bytes[0], pld_bytes[1], pld_bytes[2], pld_bytes[3]);
			fclose(fp);
			break;
		}
	}

	if (pld_size != 0 && pld_size != 0x11223344) {
		if (loaded)
			return 1;
		return 2;
	}
	return 0;
}

#define SPYD2SETUP_H
#endif /* SPYD2SETUP_H */
