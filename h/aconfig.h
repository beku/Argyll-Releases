
#ifndef __CONFIG_H__
#define __CONFIG_H__

/* General project wide configuration */


/* Version of Argyll release */
/* Bug fix = 4 bits */
/* minor number = 8 bits */
/* major number = 8 bits */

#define ARGYLL_VERSION 0x01030
#define ARGYLL_VERSION_STR "1.3.0"

/* Maximum file path length */
#define MAXNAMEL 1024

#define ENABLE_SERIAL	/* Enable access to serial instruments */
#define ENABLE_USB		/* Enable access to USB instruments using libusb */
#if defined(UNIX) && !defined(__APPLE__)
#define USE_UCMM		/* Enable the Unix micro CMM */
#endif

#endif /* __CONFIG_H__ */
