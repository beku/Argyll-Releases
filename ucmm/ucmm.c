
/* 
 * Unix micro-cmm to manage X11 display
 * calibration and profile loading.
 */

/*************************************************************************
 Copyright 2008 Graeme W. Gill

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.

 *************************************************************************/

/* We use libjcnf to store the ICC profile association with particular displays */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#ifndef NT
# include <unistd.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include "icc.h"
#include "xdg_bds.h"
#include "jcnf.h"
#include "ucmm.h"

#undef DEBUG

#define CONFIG_FILE "color.jcnf"
#define PROFILE_DIR "color/icc/devices/display"

#ifdef DEBUG
# define errout stderr
# define debug(xx)	fprintf(errout, xx )
# define debug2(xx)	fprintf xx
#else
# define debug(xx)
# define debug2(xx)
#endif

static unsigned int fnv_32_buf(void *buf, size_t len);

#ifdef NT

/* Given the path to a file, ensure that all the parent directories */
/* are created. return nz on error */
static int mkdirs(char *path) {
	struct _stat sbuf;
	char *pp = path;

	if (*pp != '\000'		/* Skip drive number */
		&& ((*pp >= 'a' && *pp <= 'z') || (*pp >= 'A' && *pp <= 'Z'))
	    && pp[1] == ':')
		pp += 2;
	if (*pp == '/')
		pp++;			/* Skip root directory */
	for (;pp != NULL && *pp != '\000';) {
		if ((pp = strchr(pp, '/')) != NULL) {
			*pp = '\000';
			if (_stat(path,&sbuf) != 0)
			{
				if (_mkdir(path) != 0)
					return 1;
			}
			*pp = '/';
			pp++;
		}
	}
	return 0;
}

#else

/* Given the path to a file, ensure that all the parent directories */
/* are created. return nz on error */
static int mkpdirs(char *path) {
	struct stat sbuf;
	char *pp = path;
	mode_t mode = 0700;		/* Default directory mode */

	if (*pp == '/')
		pp++;			/* Skip root directory */
	for (;pp != NULL && *pp != '\000';) {
		if ((pp = strchr(pp, '/')) != NULL) {
			*pp = '\000';
			if (stat(path,&sbuf) != 0) {
				if (mkdir(path, mode) != 0)
					return 1;
			} else
				mode = sbuf.st_mode;
			*pp = '/';
			pp++;
		}
	}
	return 0;
} 
#endif /* !NT */

/* Given a block of binary, convert it to upper case hexadecimal, */
/* with a 0x prefix. Free the buffer returned. */
/* Return NULL on error */
static char *buf2hex(unsigned char *buf, int len) {
	char *s;
	int i;
	
	char hex[17] = "0123456789ABCDEF";

	if ((s = malloc(len * 2 + 3)) == NULL)
		return NULL;

	s[0] = '0';
	s[1] = 'x';

	for (i = 0; i < len; i++) {
		s[2 + i * 2 + 0] = hex[(buf[i] >> 4) & 0xf];
		s[2 + i * 2 + 1] = hex[buf[i] & 0xf];
	}
	s[2 + i * 2 + 0] = '\000';

	return s;
}


/* Install a profile for a given monitor */
/* Either EDID or display_name may be NULL, but not both. */
/* Any existing association is overwritten. Installed profiles */
/* are not deleted. */
ucmm_error ucmm_install_monitor_profile(
	ucmm_scope scope,		/* Scope of instalation */
	unsigned char *edid,	/* Primary device identifier, NULL if none. */
	int edid_len,			/* Length of edid data */
	char *display_name,		/* Fall back device association, */
	                        /* the X11 display name */
	char *profile			/* Path to profile to be installed. */
) {
	char *config_file = CONFIG_FILE;
	char *profile_dir = PROFILE_DIR;
	char *conf_name = NULL;		/* Configuration path to use */
	char *data_name = NULL;		/* Data path to use */
	char *dprof = NULL;			/* Destination for profile */
	unsigned int edid_hash = 0;

	if (edid != NULL)
		edid_hash = fnv_32_buf(edid, edid_len);

	debug2((errout,"ucmm_install_monitor_profile called with profile '%s', edid 0x%x, disp '%s'\n",profile,edid_hash,display_name));

	/* Verify that we've been given a suitable ICC profile */
	/* And read it into a memory buffer */
	{
		icmFile *fp;
		icc *icco;

		if ((fp = new_icmFileStd_name(profile,"r")) == NULL) {
			debug2((errout,"Unable to ope file '%s'\n",profile));
			return ucmm_invalid_profile;
		}

		if ((icco = new_icc()) == NULL) {
			debug2((errout,"new_icc() failed\n"));
			fp->del(fp);
			return ucmm_invalid_profile;
		}

		if (icco->read(icco,fp,0) != 0) {
			debug2((errout,"icc read of '%s' failed\n",profile));
			icco->del(icco);
			fp->del(fp);
			return ucmm_invalid_profile;
		}

		if (icco->header->deviceClass != icSigDisplayClass
		 || icco->header->colorSpace != icSigRgbData) {
			debug2((errout,"profile '%s' isn't an RGB display profile\n",profile));
			icco->del(icco);
			fp->del(fp);
			return ucmm_invalid_profile;
		}
		icco->del(icco);
		fp->del(fp);
	}

	debug2((errout,"verified profile OK\n"));

	/* Locate the directories where the config file is, */
	/* and where we should copy the profile to. */
	{
		int npaths;
		xdg_error er;
		char *data_pathfile;		/* Path & name of destination */
		char **paths;
		char *tt;

		if ((npaths = xdg_bds(&er, &paths, xdg_conf, xdg_write, 
		                     scope == ucmm_local_system ? xdg_local : xdg_user,
		                     config_file)) == 0) {
			return ucmm_open_config;
		}
		if ((conf_name = strdup(paths[0])) == NULL) {
			xdg_free(paths, npaths);
			return ucmm_resource;
		}
		xdg_free(paths, npaths);
		
		/* Combined sub-path and profile name */
		if ((data_pathfile = malloc(strlen(profile_dir) + 1 + strlen(profile) + 1)) == NULL)
			return ucmm_resource;
		strcpy(data_pathfile, profile_dir);

		if (strlen(data_pathfile) > 1 && data_pathfile[strlen(data_pathfile)-1] != '/')
			strcat(data_pathfile, "/");

		if ((tt = strrchr(profile, '/')) != NULL)	/* Get base name of profile */
			tt++;
		else
			tt = profile;
		strcat(data_pathfile, tt);

		if ((npaths = xdg_bds(&er, &paths, xdg_conf, xdg_write, 
		                     scope == ucmm_local_system ? xdg_local : xdg_user,
		                     data_pathfile)) == 0) {
			free(data_pathfile);
			free(conf_name);
			return ucmm_open_config;
		}
		free(data_pathfile);
		if ((data_name = strdup(paths[0])) == NULL) {
			free(conf_name);
			xdg_free(paths, npaths);
			return ucmm_resource;
		}
		xdg_free(paths, npaths);
	}

	debug2((errout,"config file = '%s'\n",conf_name));
	debug2((errout,"data file = '%s'\n",data_name));
	
	/* Copy the profile to the destination */
	{
		FILE *fp;
		unsigned char *pdata;
		unsigned long psize;

		/* Read in the ICC profile, then set the X11 atom value */
#if defined(O_BINARY) || defined(_O_BINARY)
		if ((fp = fopen(profile,"rb")) == NULL)
#else
		if ((fp = fopen(profile,"r")) == NULL)
#endif
		{
			debug2((errout,"Can't open file '%s'\n",profile));
			free(conf_name);
			free(data_name);
			return ucmm_profile_copy;
		}

		/* Figure out how big it is */
		if (fseek(fp, 0, SEEK_END)) {
			debug2((errout,"Seek '%s' to EOF failed\n",profile));
			free(conf_name);
			free(data_name);
			return ucmm_profile_copy;
		}
		psize = (unsigned long)ftell(fp);
	
		if (fseek(fp, 0, SEEK_SET)) {
			debug2((errout,"Seek '%s' to SOF failed\n",profile));
			free(conf_name);
			free(data_name);
			return ucmm_profile_copy;
		}
	
		if ((pdata = (unsigned char *)malloc(psize)) == NULL) {
			debug2((errout,"Failed to allocate buffer for profile '%s'\n",profile));
			free(conf_name);
			free(data_name);
			return ucmm_profile_copy;
		}
	
		if (fread(pdata, 1, psize, fp) != psize) {
			debug2((errout,"Failed to read profile '%s' into buffer\n",profile));
			free(conf_name);
			free(data_name);
			return ucmm_profile_copy;
		}
		
		fclose(fp);

		/* Write the profile to its location */
		if (mkpdirs(data_name)) {
			debug2((errout,"Can't create directories for file '%s'\n",data_name));
			free(conf_name);
			free(data_name);
			return ucmm_profile_copy;
		}
#if defined(O_BINARY) || defined(_O_BINARY)
		if ((fp = fopen(data_name,"wb")) == NULL)
#else
		if ((fp = fopen(data_name,"w")) == NULL)
#endif
		{
			debug2((errout,"Can't create file '%s'\n",data_name));
			free(conf_name);
			free(data_name);
			return ucmm_profile_copy;
		}

		if (fwrite(pdata, 1, psize, fp) != psize) {
			debug2((errout,"Failed to write profile '%s' into buffer\n",data_name));
			free(conf_name);
			free(data_name);
			return ucmm_profile_copy;
		}

		if (fclose(fp) != 0) {
			debug2((errout,"Failed to close profile '%s' into buffer\n",data_name));
			free(conf_name);
			free(data_name);
			return ucmm_profile_copy;
		}
	}

	debug2((errout,"profile copied OK\n"));

	/* Update the config file */
	{
		jc_error ev;
		jcnf *jc;
		char keyn1[100];
		char keyn2[100];
		char *mname;			/* Name of key to match to */
		char *mval;				/* Value to match */
		int ix = 0;
		int recno = -1;		/* Number of the last record read */

		/* Open the configuration file for modification */
		if (mkpdirs(conf_name)) {
			debug2((errout,"Can't create directories for file '%s'\n",conf_name));
			free(conf_name);
			free(data_name);
			return ucmm_open_config;
		}

		if ((jc = new_jcnf(&ev, conf_name, jc_modify, jc_create)) == NULL) {
			debug2((errout,"new_jcnf '%s' failed with error %d\n",conf_name,ev));
			free(conf_name);
			free(data_name);
			return ucmm_open_config;
		}
	
		/* if EDID supplied, Locate a matching EDID */
		if (edid != NULL) {
			mname = "EDID";
			if ((mval = buf2hex(edid, edid_len)) == NULL) {
				jc->del(jc);
				free(conf_name);
				free(data_name);
				return ucmm_resource;
			}

		/* Else fall back to X11 display name and screen */
		} else {
			if (display_name == NULL) {
				jc->del(jc);
				free(conf_name);
				free(data_name);
				return ucmm_no_edid_or_display;
			}
			mname = "NAME";
			if ((mval = strdup(display_name)) == NULL) {
				jc->del(jc);
				free(conf_name);
				free(data_name);
				return ucmm_resource;
			}
		}

		debug2((errout,"Searching for %s = '%s'\n",mname,mval));
		for (;;ix++) {
			char *key, *pp;
			jc_type type;
			unsigned char *data;
			size_t dataSize;
			int ii;

			if ((ev = jc->locate_key(jc, &ix, "devices/display/", 0, 0)) != jc_ok
			 || (ev = jc->get_key(jc, ix, &key, &type, &data, &dataSize, NULL)) != jc_ok) {
				if (ev == jc_ix_oorange) {
					break;
				}
				debug2((errout,"jcnf locate/get_key failed with error %d\n",ev));
				free(mval);
				jc->del(jc);
				free(conf_name);
				free(data_name);
				return ucmm_open_config;
			}
			
			if ((pp = jc_get_nth_elem(key, 2)) == NULL) {
				continue;
			}
			if ((ii = atoi(pp)) == 0) {
				free(pp);
				continue;
			}
			if (ii > recno)		/* Track biggest, so we know what to create next */
				recno = ii;
			if ((pp = jc_get_nth_elem(key, 3)) != NULL && strcmp(pp, mname) == 0 && type == jc_string && strcmp(data, mval) == 0) {
				/* Found matching record */
				free(pp);
				break;
			}
			if (pp != NULL)
				free(pp);
		}
		
		/* Create a new record */
		if (ev == jc_ix_oorange) {
			recno++;				/* Make it the next index */
			if (recno <= 0)
				recno = 1;
			debug2((errout, "Adding a new record %d\n",recno));
		} else {
			debug2((errout, "Replacing record %d\n",recno));
		}

		/* Write the record */
		sprintf(keyn1, "devices/display/%d/%s", recno, mname);
		sprintf(keyn2, "devices/display/%d/ICC_PROFILE", recno);
		if ((ev = jc->set_key(jc, -1, keyn1, jc_string, mval, strlen(mval)+1, NULL)) != jc_ok
		 || (ev = jc->set_key(jc, -1, keyn2, jc_string, data_name, strlen(data_name)+1, NULL)) != jc_ok) {
			debug2((errout,"jcnf set_key failed with error %d\n",ev));
			free(mval);
			jc->del(jc);
			free(conf_name);
			free(data_name);
			return ucmm_set_config;
		}
		free(mval);

		/* write to record the EDID or display name and the profile path */
		if ((ev = jc->update(jc)) != 0) {
			debug2((errout,"jcnf write to '%s' failed with error %d\n",conf_name,ev));
			jc->del(jc);
			free(conf_name);
			free(data_name);
			return ucmm_save_config;
		}
		debug2((errout,"Updated config file '%s'\n",conf_name));

		/* We're done with this */
		jc->del(jc);
		free(conf_name);
		free(data_name);
	}
	debug2((errout,"ucmm done profile install\n"));
	return ucmm_ok;
}

/* Un-install a profile for a given monitor. */
/* Either EDID or display_name may be NULL, but not both. */
/* The monitor is left with no profile association. If a profile */
/* name is provided and matches the one that was associated with */
/* the monitor, and has no other association, then it will be deleted */
/* from the data directory. */
/* Return an error code */
ucmm_error ucmm_uninstall_monitor_profile(
	ucmm_scope scope,		/* Scope of instalation */
	unsigned char *edid,	/* Primary device identifier, NULL if none. */
	int edid_len,			/* Length of edid data */
	char *display_name,		/* Fall back device association, */
	char *profile			/* Base name of profile to be deleted. NULL if not to be deleted. */
) {
	char *config_file = CONFIG_FILE;
	char *profile_dir = PROFILE_DIR;
	char *conf_name = NULL;		/* Configuration path to use */
	char *data_name = NULL;		/* Data path to use */
	char *dprof = NULL;			/* Destination for profile */
	unsigned int edid_hash = 0;

	if (edid != NULL)
		edid_hash = fnv_32_buf(edid, edid_len);

	debug2((errout,"ucmm_uninstall_monitor_profile called with profile '%s', edid 0x%x, disp '%s'\n",profile,edid_hash,display_name));

	/* Locate the directories where the config file is, */
	/* and where the profile should be too. */
	{
		int npaths;
		xdg_error er;
		char *data_pathfile;		/* Path & name of destination */
		char **paths;
		char *tt;

		if ((npaths = xdg_bds(&er, &paths, xdg_conf, xdg_read, 
		                     scope == ucmm_local_system ? xdg_local : xdg_user,
		                     config_file)) == 0) {
			return ucmm_open_config;
		}
		if ((conf_name = strdup(paths[0])) == NULL) {
			xdg_free(paths, npaths);
			return ucmm_resource;
		}
		xdg_free(paths, npaths);
		
		if (profile != NULL) {
			/* Combined sub-path and profile name */
			if ((data_pathfile = malloc(strlen(profile_dir) + 1 + strlen(profile) + 1)) == NULL)
				return ucmm_resource;
			strcpy(data_pathfile, profile_dir);

			if (strlen(data_pathfile) > 1 && data_pathfile[strlen(data_pathfile)-1] != '/')
				strcat(data_pathfile, "/");

			if ((tt = strrchr(profile, '/')) != NULL)	/* Get base name of profile */
				tt++;
			else
				tt = profile;
			strcat(data_pathfile, tt);
	
			if ((npaths = xdg_bds(&er, &paths, xdg_conf, xdg_read, 
			                     scope == ucmm_local_system ? xdg_local : xdg_user,
			                     data_pathfile)) == 0) {
				free(data_pathfile);
				free(conf_name);
				return ucmm_open_config;
			}
			free(data_pathfile);
			if ((data_name = strdup(paths[0])) == NULL) {
				free(conf_name);
				xdg_free(paths, npaths);
				return ucmm_resource;
			}
			xdg_free(paths, npaths);
		}
	}

	debug2((errout,"config file = '%s'\n",conf_name));
	if (data_name != NULL)
		debug2((errout,"data file = '%s'\n",data_name));
	
	/* Get the config file */
	{
		jc_error ev;
		jcnf *jc;
		char keyn1[100];
		char *mname;			/* Name of key to match to */
		char *mval;				/* Value to match */
		int ix;
		int recno = -1;		/* Number of the last record read */

		/* Open the configuration file for modification */
		if (mkpdirs(conf_name)) {
			debug2((errout,"Can't create directories for file '%s'\n",conf_name));
			free(conf_name);
			if (data_name != NULL)
				free(data_name);
			return ucmm_open_config;
		}

		if ((jc = new_jcnf(&ev, conf_name, jc_modify, jc_create)) == NULL) {
			debug2((errout,"new_jcnf '%s' failed with error %d\n",conf_name,ev));
			free(conf_name);
			if (data_name != NULL)
				free(data_name);
			return ucmm_open_config;
		}
	
		/* if EDID supplied, Locate a matching EDID */
		if (edid != NULL) {
			mname = "EDID";
			if ((mval = buf2hex(edid, edid_len)) == NULL) {
				jc->del(jc);
				free(conf_name);
				if (data_name != NULL)
					free(data_name);
				return ucmm_resource;
			}

		/* Else fall back to X11 display name and screen */
		} else {
			if (display_name == NULL) {
				jc->del(jc);
				free(conf_name);
				if (data_name != NULL)
					free(data_name);
				return ucmm_no_edid_or_display;
			}
			mname = "NAME";
			if ((mval = strdup(display_name)) == NULL) {
				jc->del(jc);
				free(conf_name);
				if (data_name != NULL)
					free(data_name);
				return ucmm_resource;
			}
		}

		debug2((errout,"Searching for %s = '%s'\n",mname,mval));
		for (ix = 0;;ix++) {
			char *key, *pp;
			jc_type type;
			unsigned char *data;
			size_t dataSize;
			int ii;

			if ((ev = jc->locate_key(jc, &ix, "devices/display/", 0, 0)) != jc_ok
			 || (ev = jc->get_key(jc, ix, &key, &type, &data, &dataSize, NULL)) != jc_ok) {
				if (ev == jc_ix_oorange) {
					break;
				}
				debug2((errout,"jcnf locate/get_key failed with error %d\n",ev));
				free(mval);
				jc->del(jc);
				free(conf_name);
				if (data_name != NULL)
					free(data_name);
				return ucmm_open_config;
			}
			
			if ((pp = jc_get_nth_elem(key, 2)) == NULL) {
				continue;
			}
			if ((ii = atoi(pp)) == 0) {
				free(pp);
				continue;
			}
			if (ii > recno)		/* Track biggest, so we know what to create next */
				recno = ii;
			if ((pp = jc_get_nth_elem(key, 3)) != NULL && strcmp(pp, mname) == 0 && type == jc_string && strcmp(data, mval) == 0) {
				/* Found matching record */
				free(pp);
				break;
			}
			if (pp != NULL)
				free(pp);
		}
		
		if (ev == jc_ix_oorange) {
			debug2((errout,"No matching display was found\n"));
			free(mval);
			jc->del(jc);
			free(conf_name);
			if (data_name != NULL)
				free(data_name);
			return ucmm_monitor_not_found;
			/* (Should we delete the file anyway ???) */
		}
		free(mval);

		debug2((errout,"Deleting record %d key '%s'\n",recno,keyn1));

		/* Delete the record */
		sprintf(keyn1, "devices/display/%d/", recno);

		for (ix = -1;;ix--) {
			if ((ev = jc->locate_key(jc, &ix, keyn1, 0, 1)) == jc_ok) {
				if ((ev = jc->delete_key(jc, ix, NULL)) != jc_ok) {
					debug2((errout,"jcnf delete_key failed with error %d\n",ev));
					jc->del(jc);
					free(conf_name);
					if (data_name != NULL)
						free(data_name);
					return ucmm_delete_key;
				}
			} else {
				if (ev == jc_ix_oorange) {
					break;
				}
				debug2((errout,"jcnf locate/get_key failed with error %d\n",ev));
				jc->del(jc);
				free(conf_name);
				if (data_name != NULL)
					free(data_name);
				return ucmm_open_config;
			}
		}

		if (data_name != NULL) {
			/* See if the profile is used by any other device */

			debug2((errout, "Searching for any reference to profile '%s'\n",data_name));
			for (ix = 0;;ix++) {
				char *key, *pp;
				jc_type type;
				unsigned char *data;
				size_t dataSize;

				if ((ev = jc->locate_key(jc, &ix, "devices/display/", 0, 0)) != jc_ok
				 || (ev = jc->get_key(jc, ix, &key, &type, &data, &dataSize, NULL)) != jc_ok) {
					if (ev == jc_ix_oorange) {
						break;
					}
					debug2((errout,"jcnf locate/get_key failed with error %d\n",ev));
					jc->del(jc);
					free(conf_name);
					if (data_name != NULL)
						free(data_name);
					return ucmm_access_config;
				}
				if ((pp = jc_get_nth_elem(key, 3)) == NULL)
					continue;
				if (strcmp(pp,"ICC_PROFILE") != 0
				 || type != jc_string
				 || strcmp(data, data_name) != 0) {
					free(pp);
					continue;
				}
				free(pp);
				break;
			}
			/* If not, delete the file */
			if (ev == jc_ix_oorange) {
				debug2((errout,"Deleting profile '%s'\n",data_name));
				if (unlink(data_name) != 0) {
					debug2((errout,"ucmm unlink '%s' failed\n",data_name));
					jc->del(jc);
					free(conf_name);
					if (data_name != NULL)
						free(data_name);
					return ucmm_access_config;
				}
			}
		}

		/* Update the config */
		if ((ev = jc->update(jc)) != 0) {
			debug2((errout,"jcnf write to '%s' failed with error %d\n",conf_name,ev));
			jc->del(jc);
			free(conf_name);
			if (data_name != NULL)
				free(data_name);
			return ucmm_save_config;
		}
		debug2((errout,"Updated config file '%s'\n",conf_name));

		/* We're done with this */
		jc->del(jc);
		free(conf_name);
		if (data_name != NULL)
			free(data_name);
	}
	debug2((errout,"ucmm done profile un-install\n"));
	return ucmm_ok;
}

/* Get an associated monitor profile. */
/* Return ucmm_no_profile if there is no installed profile for this */
/* monitor. */
/* Return an error code */
ucmm_error ucmm_get_monitor_profile(
	unsigned char *edid,	/* Primary device identifier, NULL if none. */
	int edid_len,			/* Length of edid data */
	char *display_name,	    /* Fall back device association, */
	char **profile		    /* Return path to profile. free() afterwards. */
) {
	int scope;
	char *config_file = "color.jcnf";
	char *conf_name = NULL;		/* Configuration path to use */
	unsigned int edid_hash = 0;

	if (edid != NULL)
		edid_hash = fnv_32_buf(edid, edid_len);

	debug2((errout,"ucmm_get_monitor_profile called edid 0x%x, disp '%s'\n",edid_hash,display_name));

	/* Look at user then local system scope */
	for (scope = 0; scope < 2; scope++) {
 
		/* Locate the directories where the config file is, */
		{
			int npaths;
			xdg_error er;
			char **paths;
			char *tt;

			if ((npaths = xdg_bds(&er, &paths, xdg_conf, xdg_read, 
			                     scope == ucmm_local_system ? xdg_local : xdg_user,
			                     config_file)) == 0) {
				continue;
			}
			if ((conf_name = strdup(paths[0])) == NULL) {
				xdg_free(paths, npaths);
				return ucmm_resource;
			}
			xdg_free(paths, npaths);
		}

		/* Get the config file */
		{
			jc_error ev;
			jcnf *jc;
			char keyn1[100];
			char *mname;			/* Name of key to match to */
			char *mval;				/* Value to match */
			int ix;
			int recno = -1;		/* Number of the last record read */
			char *key, *pp;
			jc_type type;
			unsigned char *data;
			size_t dataSize;

			/* Open the configuration file for reading */
			if ((jc = new_jcnf(&ev, conf_name, jc_read, jc_no_create)) == NULL) {
				debug2((errout,"new_jcnf '%s' failed with error %d\n",conf_name,ev));
				continue;		/* Try the next scope */
			}
		
			/* if EDID supplied, Locate a matching EDID */
			if (edid != NULL) {
				mname = "EDID";
				if ((mval = buf2hex(edid, edid_len)) == NULL) {
					debug2((errout,"buf2jex  failed\n"));
					jc->del(jc);
					free(conf_name);
					return ucmm_resource;
				}

			/* Else fall back to X11 display name and screen */
			} else {
				if (display_name == NULL) {
					debug2((errout,"No EDID and display name  failed\n"));
					jc->del(jc);
					free(conf_name);
					return ucmm_no_edid_or_display;
				}
				mname = "NAME";
				if ((mval = strdup(display_name)) == NULL) {
					debug2((errout,"strdup failed\n"));
					jc->del(jc);
					free(conf_name);
					return ucmm_resource;
				}
			}

			debug2((errout,"Searching for %s = '%s'\n",mname,mval));
			for (ix = 0;;ix++) {
				int ii;

				if ((ev = jc->locate_key(jc, &ix, "devices/display/", 0, 0)) != jc_ok
				 || (ev = jc->get_key(jc, ix, &key, &type, &data, &dataSize, NULL)) != jc_ok) {
					if (ev == jc_ix_oorange) {
						break;
					}
					debug2((errout,"jcnf locate/get_key failed with error %d\n",ev));
					free(mval);
					jc->del(jc);
					free(conf_name);
					return ucmm_open_config;
				}
				
				if ((pp = jc_get_nth_elem(key, 2)) == NULL) {
					continue;
				}
				if ((ii = atoi(pp)) == 0) {
					free(pp);
					continue;
				}
				if (ii > recno)		/* Track biggest, so we know what to create next */
					recno = ii;
				if ((pp = jc_get_nth_elem(key, 3)) != NULL && strcmp(pp, mname) == 0 && type == jc_string
				  && strcmp(data, mval) == 0) {
					/* Found matching record */
					free(pp);
					break;
				}
				if (pp != NULL)
					free(pp);
			}
			
			if (ev == jc_ix_oorange) {
				debug2((errout,"No matching display was found\n"));
				continue;	/* On to the next scope */
			}
			free(mval);

			/* Get the profile path from the  record */
			sprintf(keyn1, "devices/display/%d/ICC_PROFILE", recno);
			key = keyn1;
			debug2((errout,"Looking up record %d key '%s'\n",recno,keyn1));

			if ((ev = jc->get_key(jc, -1, &key, &type, &data, &dataSize, NULL)) != jc_ok
			 || type != jc_string) {
				debug2((errout,"jcnf locate/get_key failed with error %d\n",ev));
				jc->del(jc);
				free(conf_name);
				if (ev == jc_ix_oorange) {
					continue;			/* try the next config */
				}
				return ucmm_access_config;
			}
			if ((*profile = strdup(data)) == NULL) {
				debug2((errout,"jcnf get_key malloc failed\n"));
				jc->del(jc);
				free(conf_name);
				return ucmm_resource;
			}
				
			/* We're done with this */
			jc->del(jc);
			free(conf_name);
			return ucmm_ok;
			debug2((errout,"Returning current profile '%s'\n",data));
		}
	}
	debug2((errout,"Failed to find a current profile\n"));

	return ucmm_no_profile; 
}
	

/* Return an ASCII error message string interpretation of an error number */
char *ucmm_error_string(ucmm_error erno) { 
	
	switch (erno) {
		case ucmm_ok:
			return "OK";
		case ucmm_resource:
			return "Resource failure (e.g. out of memory)";
		case ucmm_invalid_profile:
			return "Profile is not a valid display ICC profile";
		case ucmm_no_profile:
			return "There is no associated profile";
		case ucmm_no_home:
			return "There is no HOME environment variable defined";
		case ucmm_no_edid_or_display:
			return "There is no edid or display name";
		case ucmm_profile_copy:
			return "There was an error copying the profile";
		case ucmm_open_config:
			return "There was an error opening the config file";
		case ucmm_access_config:
			return "There was an error accessing the config information";
		case ucmm_set_config:
			return "There was an error setting the config file";
		case ucmm_save_config:
			return "There was an error saving the config file";
		case ucmm_monitor_not_found:
			return "The EDID or display wasn't matched";
		case ucmm_delete_key:
			return "Delete_key failed";
		case ucmm_delete_profile:
			return "Delete_key failed";
	}
	return "Unknown error number";
}


/* ============================================================= */

/*
 * hash_32 - 32 bit Fowler/Noll/Vo hash code
 *
 * @(#) $Revision: 1.8 $
 * @(#) $Id: hash_32.c,v 1.8 2003/10/03 20:38:13 chongo Exp $
 * @(#) $Source: /usr/local/src/cmd/fnv/RCS/hash_32.c,v $
 *
 ***
 *
 * Fowler/Noll/Vo hash
 *
 * The basis of this hash algorithm was taken from an idea sent
 * as reviewer comments to the IEEE POSIX P1003.2 committee by:
 *
 *	  Phong Vo (http://www.research.att.com/info/kpv/)
 *	  Glenn Fowler (http://www.research.att.com/~gsf/)
 *
 * In a subsequent ballot round:
 *
 *	  Landon Curt Noll (http://www.isthe.com/chongo/)
 *
 * improved on their algorithm.  Some people tried this hash
 * and found that it worked rather well.  In an EMail message
 * to Landon, they named it the ``Fowler/Noll/Vo'' or FNV hash.
 *
 * FNV hashes are designed to be fast while maintaining a low
 * collision rate. The FNV speed allows one to quickly hash lots
 * of data while maintaining a reasonable collision rate.  See:
 *
 *	  http://www.isthe.com/chongo/tech/comp/fnv/index.html
 *
 * for more details as well as other forms of the FNV hash.
 ***
 *
 * NOTE: The FNV-0 historic hash is not recommended.  One should use
 *	 the FNV-1 hash instead.
 *
 * To use the 32 bit FNV-0 historic hash, pass FNV0_32_INIT as the
 * unsigned int hashval argument to fnv_32_buf() or fnv_32_str().
 *
 * To use the recommended 32 bit FNV-1 hash, pass FNV1_32_INIT as the
 * unsigned int hashval argument to fnv_32_buf() or fnv_32_str().
 *
 ***
 *
 * Please do not copyright this code.  This code is in the public domain.
 *
 * LANDON CURT NOLL DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
 * INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO
 * EVENT SHALL LANDON CURT NOLL BE LIABLE FOR ANY SPECIAL, INDIRECT OR
 * CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
 * USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 * PERFORMANCE OF THIS SOFTWARE.
 *
 * By:
 *	chongo <Landon Curt Noll> /\oo/\
 *	  http://www.isthe.com/chongo/
 *
 * Share and Enjoy!	:-)
 */

/*
 * 32 bit magic FNV-0 and FNV-1 prime
 */
#define FNV_32_PRIME ((unsigned int)0x01000193)

#define FNV1_32_INIT ((unsigned int)0x811c9dc5)

/*
 * fnv_32_buf - perform a 32 bit Fowler/Noll/Vo hash on a buffer
 *
 * input:
 *	buf	- start of buffer to hash
 *	len	- length of buffer in octets
 *	hval	- previous hash value or 0 if first call
 *
 * returns:
 *	32 bit hash as a static hash type
 *
 * NOTE: To use the recommended 32 bit FNV-1 hash, use FNV1_32_INIT as the hval
 *	 argument on the first call to either fnv_32_buf() or fnv_32_str().
 */
static unsigned int
fnv_32_buf_cont(void *buf, size_t len, unsigned int hval)
{
	unsigned char *bp = (unsigned char *)buf;	/* start of buffer */
	unsigned char *be = bp + len;		/* beyond end of buffer */

	/*
	 * FNV-1 hash each octet in the buffer
	 */
	while (bp < be) {

	/* multiply by the 32 bit FNV magic prime mod 2^32 */
#if defined(NO_FNV_GCC_OPTIMIZATION)
	hval *= FNV_32_PRIME;
#else
	hval += (hval<<1) + (hval<<4) + (hval<<7) + (hval<<8) + (hval<<24);
#endif

	/* xor the bottom with the current octet */
	hval ^= (unsigned int)*bp++;
	}

	/* return our new hash value */
	return hval;
}

static unsigned int
fnv_32_buf(void *buf, size_t len) {
	return fnv_32_buf_cont(buf, len, FNV1_32_INIT);
}

