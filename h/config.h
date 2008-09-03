
#ifndef __CONFIG_H__
#define __CONFIG_H__

/* General project wide configuration */


/* Version of Argyll release */
/* Bug fix = 4 bits */
/* minor number = 8 bits */
/* major number = 8 bits */

#define ARGYLL_VERSION 0x01003
#define ARGYLL_VERSION_STR "1.0.3"

/* Maximum file path length */
#define MAXNAMEL 512

/* Enable access to USB instruments using libusb */
#define ENABLE_USB

#endif /* __CONFIG_H__ */
