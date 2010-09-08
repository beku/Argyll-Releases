
 /* XDG Base Directory Specifications support library. */
 /* Implements equivalent cross platform functionality too. */

/* 
 * Argyll Color Correction System
 *
 * Author: Graeme W. Gill
 * Date:   28/7/2010
 *
 * Copyright 2010 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
	This function provides support for the XDG Base Directory Specifications
	in a cross platform compatible way.

	[ Note that for MSWin each path in a set is separated by a ';' character.
	  and that DATA and CONF will be in the same directory. ]

	The following paths are use for each of the 5 XDG concepts, listed in order
	of priority:

		Per user application related data.

		Per user application configuration settings.

		Per user application cache storage area.

		Local system wide application related data.

		Local system wide application configuration settings.

	Unix:
		$XDG_DATA_HOME
		$HOME/.local/share

		$XDG_CONF_HOME
		$HOME/.config

		$XDG_CACHE_HOME
		$HOME/.cache

		$XDG_DATA_DIRS
		/usr/local/share:/usr/share

		$XDG_CONF_DIRS
		/etc/xdg

	OS X:
		$XDG_DATA_HOME
		$HOME/Library

		$XDG_CONF_HOME
		$HOME/Library/Preferences

		$XDG_CACHE_HOME
		$HOME/Library/Caches

		$XDG_DATA_DIRS
		/Library

		$XDG_CONF_DIRS
		/Library/Preferences

	MSWin:
		$XDG_DATA_HOME
		$APPDATA
		$HOME/.local/share
		$USERPROFILE/Application Data

		$XDG_CONF_HOME
		$APPDATA
		$HOME/.config
		$USERPROFILE/Application Data

		$XDG_CACHE_HOME
		$APPDATA/Cache
		$HOME/.cache
		$USERPROFILE/Application Data/Cache

		$XDG_DATA_DIRS
		$ALLUSERSPROFILE

		$XDG_CONF_DIRS
		$ALLUSERSPROFILE
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#ifndef NT
#include <unistd.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include "xdg_bds.h"

#undef DEBUG
#define dbgo stderr

#ifdef DEBUG
#define DBG(xxx) fprintf xxx ;
#else
#define DBG(xxx) 
#endif	/* DEBUG */

#ifdef NT
#define SSEP ';'		/* Since ':' is used for drive letter */
#define SSEPS ";"
#else
#define SSEP ':'
#define SSEPS ":"
#endif

#ifdef NT
# define stat _stat
# define mode_t int
# define mkdir(A,B) _mkdir(A)
# define mputenv _putenv
# define unlink _unlink
# define rmdir _rmdir
#else
/* UNIX putenv is a pain.. */
static void mputenv(char *ss) {
	int ll = strlen(ss);
	ss = strdup(ss);
	if (ll > 0 && ss[ll-1]== '=') {
		ss[ll-1] = '\000';
		unsetenv(ss); 
	} else {
		putenv(ss);
	}
}
#endif

/* Allocate a copy of the string, and normalize the */
/* path separator to '/' */

/* Append a string. Free in. Return NULL on error. */
static char *append(char *in, char *app) {
	char *rv;

	if ((rv = malloc(strlen(in) + strlen(app) + 1)) == NULL) {
		free(in);
		return NULL;
	}
	strcpy(rv, in);
	strcat(rv, app);
	free(in);

	return rv;
}

/* Append a ':'/';' then a string. Free in. Return NULL on error. */
static char *cappend(char *in, char *app) {
	int inlen;
	char *rv;

	inlen = strlen(in);

	if ((rv = malloc(inlen + 1 + strlen(app) + 1)) == NULL) {
		free(in);
		return NULL;
	}
	strcpy(rv, in);
	if (inlen > 1)
		strcat(rv, SSEPS);
	strcat(rv, app);
	free(in);

	return rv;
}

/* Append a '/' then a string. Free in. Return NULL on error. */
static char *dappend(char *in, char *app) {
	int inlen;
	char *rv;

	inlen = strlen(in);

	if ((rv = malloc(inlen + 1 + strlen(app) + 1)) == NULL) {
		free(in);
		return NULL;
	}
	strcpy(rv, in);
	if (inlen > 1 && in[inlen-1] != '/')
		strcat(rv, "/");
	strcat(rv, app);
	free(in);

	return rv;
}

/* Return the full path to the given subpath for the type of storage */
/* and access required. Return NULL if there is an error. */
/* The string should be free()'s after use. */
/* XDG environment variables and the subpath are assumed to be using */
/* the '/' path separator. */
/* When "xdg_write", the necessary path to the file will be created. */
/* If we're running as sudo and are creating a user dir/file, */
/* we drop to using the underlying SUDO_UID/GID. If we are creating a */
/* local system dir/file as sudo and have dropped to the SUDO_UID/GID, */
/* then revert back to root uid/gid. */
char *xdg_bds(
	xdg_error *er,			/* Return an error code */
	xdg_storage_type st,	/* Specify the storage type */
	xdg_op_type op,			/* Operation type */
	xdg_scope sc,			/* Scope if write */
	char *pfname			/* Sub-path and file name */
) {
	char *path = NULL;

	DBG((dbgo,"xdg_bds called with st %s, op %s, sc %s, pfname '%s'\n",
	st == xdg_data ? "data" : st == xdg_conf ? "config" : st == xdg_cache ? "cache" : "unknown", 
	op == xdg_write ? "write" : op == xdg_read ? "read" : "unknown", 
	sc == xdg_user ? "user" : sc == xdg_local ? "local" : "unknown", 
	pfname))

	/* Initial, empty path */
	if ((path = strdup("")) == NULL) {
		if (er != NULL) *er = xdg_alloc;
		DBG((dbgo,"malloc error\n"))
		return NULL;
	}

	/* Create a set of ':'/';' separated search paths */

	/* User scope */
	if (op == xdg_read || sc == xdg_user) {
		if (st == xdg_data) {
			char *xdg, *home;
			if ((xdg = getenv("XDG_DATA_HOME")) != NULL) {
				if ((path = cappend(path, xdg)) == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
#ifdef NT
			} else if (getenv("HOME") == NULL && (xdg = getenv("APPDATA")) != NULL) {
				if ((path = cappend(path, xdg)) == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
#endif
			} else {
				if ((home = getenv("HOME")) == NULL
#ifdef NT
				  && (home = getenv("USERPROFILE")) == NULL
#endif
				) {
					if (er != NULL) *er = xdg_nohome;
					free(path);
					DBG((dbgo,"no $HOME\n"))
					return NULL;
				}
				if ((path = cappend(path, home)) == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
#ifdef NT
				if (getenv("HOME") != NULL)
					path = dappend(path, ".local/share");
				else
					path = dappend(path, "Application Data");
#else
#ifdef __APPLE__
				path = dappend(path, "Library");
#else	/* Unix, Default */
				path = dappend(path, ".local/share");
#endif
#endif
				if (path == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
			}
		} else if (st == xdg_conf) {
			char *xdg, *home;
			if ((xdg = getenv("XDG_CONF_HOME")) != NULL) {
				if ((path = cappend(path, xdg)) == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
#ifdef NT
			} else if (getenv("HOME") == NULL && (xdg = getenv("APPDATA")) != NULL) {
				if ((path = cappend(path, xdg)) == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
#endif
			} else {
				if ((home = getenv("HOME")) == NULL
#ifdef NT
				  && (home = getenv("USERPROFILE")) == NULL
#endif
				) {
					if (er != NULL) *er = xdg_nohome;
					free(path);
					DBG((dbgo,"no $HOME\n"))
					return NULL;
				}
				if ((path = cappend(path, home)) == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
#ifdef NT
				if (getenv("HOME") != NULL)
					path = dappend(path, ".config");
				else
					path = dappend(path, "Application Data");

#else
#ifdef __APPLE__
				path = dappend(path, "Library/Preferences");
#else	/* Unix, Default */
				path = dappend(path, ".config");
#endif
#endif
				if (path == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
			}
		} else if (st == xdg_cache) {
			char *xdg, *home;
			if ((xdg = getenv("XDG_CACHE_HOME")) != NULL) {
				if ((path = cappend(path, xdg)) == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
#ifdef NT
			} else if (getenv("HOME") == NULL && (xdg = getenv("APPDATA")) != NULL) {
				if ((path = cappend(path, xdg)) == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
				if ((path = dappend(path, "Cache")) == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
#endif
			} else {
				if ((home = getenv("HOME")) == NULL
#ifdef NT
				  && (home = getenv("USERPROFILE")) == NULL
#endif
				) {
					if (er != NULL) *er = xdg_nohome;
					free(path);
					DBG((dbgo,"no $HOME\n"))
					return NULL;
				}
				if ((path = cappend(path, home)) == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
#ifdef NT
				if (getenv("HOME") != NULL)
					path = dappend(path, ".cache");
				else
					path = dappend(path, "Application Data/Cache");
#else
#ifdef __APPLE__
				path = dappend(path, "Library/Caches");
#else	/* Unix, Default */
				path = dappend(path, ".cache");
#endif
#endif
				if (path == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
			}
		}
	}
	/* Local system scope */
	if (op == xdg_read || sc == xdg_local) {
		char *xdg;
		if (st == xdg_data) {
			if ((xdg = getenv("XDG_DATA_DIRS")) != NULL) {
				if ((path = cappend(path, xdg)) == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
			} else {
#ifdef NT
				/*
					QT uses $COMMON_APPDATA expected to be
						C:\Documents and Settings\All Users\Application Data\
					while others use $CommonAppData.
					Both seem poorly supported, 
				 */
				char *home;
				if ((home = getenv("ALLUSERSPROFILE")) == NULL
				) {
					if (er != NULL) *er = xdg_noalluserprofile;
					free(path);
					DBG((dbgo,"no $ALLUSERSPROFILE\n"))
					return NULL;
				}
				path = cappend(path, home);
#else
#ifdef __APPLE__
				path = cappend(path, "/Library");
#else
				path = cappend(path, "/usr/local/share:/usr/share");
#endif
#endif
				if (path == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
			}
		} else if (st == xdg_conf) {
			if ((xdg = getenv("XDG_CONF_DIRS")) != NULL) {
				if ((path = cappend(path, xdg)) == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
			} else {
#ifdef NT
				char *home;
				if ((home = getenv("ALLUSERSPROFILE")) == NULL
				) {
					if (er != NULL) *er = xdg_noalluserprofile;
					free(path);
					DBG((dbgo,"no $ALLUSERSPROFILE\n"))
					return NULL;
				}
				path = cappend(path, home);
#else
#ifdef __APPLE__
				path = cappend(path, "/Library/Preferences");
#else
				path = cappend(path, "/etc/xdg");
#endif
#endif
				if (path == NULL) {
					if (er != NULL) *er = xdg_alloc;
					DBG((dbgo,"malloc error\n"))
					return NULL;
				}
			}
		}
	}

#ifdef NT
	{
		char *cp;
		for (cp = path; *cp != '\000'; cp++) {
			if (*cp == '\\')	
				*cp = '/';
		}
	}
#endif

	DBG((dbgo,"Paths to search '%s'\n",path));

	/* Hmm. */
	if (strlen(path) == 0) {
		free(path);
		if (er != NULL) *er = xdg_nopath;
		return NULL;
	}

	{
		char *spath = NULL;
		char *cp, *ep, tc;

		/* For each search path */
		for (cp = path; *cp != '\000';) {
			char *pp;
			struct stat sbuf;
			mode_t mode;

			/* Copy search path */
			if ((ep = strchr(cp, SSEP)) == NULL)
				ep = cp + strlen(cp);
			if ((ep - cp) == 0) {
				free(path);
				if (er != NULL) *er = xdg_mallformed;
				return NULL;
			}
			if ((spath = (char *)malloc(ep - cp + 1)) == NULL) {
				free(path);
				if (er != NULL) *er = xdg_alloc;
				return NULL;
			}
			memcpy(spath, cp, ep - cp);
			spath[ep - cp] = '\000';

			/* append subpath & filename */
			if ((spath = dappend(spath, pfname)) == NULL) {
				free(path);
				if (er != NULL) *er = xdg_alloc;
				return NULL;
			} 
			DBG((dbgo,"Full path to check '%s'\n",spath));

			if (op == xdg_read) {
				struct stat sbuf;

				DBG((dbgo,"Checking path '%s'\n",spath))
				if (stat(spath,&sbuf) != 0) {
					DBG((dbgo,"Giving up on this one\n"))
					goto next_spath;
				}

			} else {	/* op == xdg_write */
				char *pp = spath;
				struct stat sbuf;
				mode_t mode = 0700;	/* Default directory mode */

				if (sc == xdg_user)
					mode = 0700;	/* Default directory mode for user */
				else
					mode = 0755;	/* Default directory mode local system shared */
#ifndef NT
				/* If we're creating a user dir/file and running as root sudo */
				if (sc == xdg_user && geteuid() == 0) {
					char *uids, *gids;
					int uid, gid;
					DBG((dbgo,"We're setting a user dir/file running as root\n"))
						
					if ((uids = getenv("SUDO_UID")) != NULL
					 && (gids = getenv("SUDO_GID")) != NULL) {
						uid = atoi(uids);
						gid = atoi(gids);
						if (setegid(gid) || seteuid(uid)) {
							DBG((dbgo,"seteuid or setegid failed\n"))
						} else {
							DBG((dbgo,"Set euid %d and egid %d\n",uid,gid))
						}
					}
				/* If setting local system dir/file and not effective root, but sudo */
				} else if (sc == xdg_local && getuid() == 0 && geteuid() != 0) {
					if (getenv("SUDO_UID") != NULL
					 && getenv("SUDO_GID") != NULL) {
						DBG((dbgo,"We're setting a local system dir/file with uid = 0 && euid != 0\n"))
						setegid(getgid());
						seteuid(getuid());
						DBG((dbgo,"Set euid %d, egid %d\n",geteuid(),getegid()))
					}
				}
#endif	/* !NT */

#ifdef NT
				if (*pp != '\000'		/* Skip drive number */
					&& ((*pp >= 'a' && *pp <= 'z') || (*pp >= 'A' && *pp <= 'Z'))
				    && pp[1] == ':')
					pp += 2;
#endif
				if (*pp == '/')
					pp++;			/* Skip root directory */

				for (;pp != NULL && *pp != '\000';) {
					if ((pp = strchr(pp, '/')) != NULL) {
						*pp = '\000';
						DBG((dbgo,"Checking path '%s'\n",spath))
						if (stat(spath,&sbuf) != 0) {
							/* Doesn't exist */
							DBG((dbgo,"Path '%s' doesn't exist\n",spath))
							if (mkdir(spath, mode) != 0) {
								DBG((dbgo,"mkdir failed - giving up on this one\n"))
								goto next_spath;
							}
						} else {
							mode = sbuf.st_mode;
						}
						*pp = '/';
						pp++;
					}
				}
			}

			/* We've found a path to return */
			if (er != NULL) *er = xdg_ok;
			free(path);
			DBG((dbgo,"Retuning '%s'\n",spath))
			return spath; 

			/* Moveto next search path */
		next_spath:;
			free(spath); spath = NULL;
			if (*ep == SSEP)
				cp = ep+1;
			else
				cp = ep;
		}
	}

	/* We've failed */
	free(path);
	if (er != NULL) *er = xdg_nopath;

	return NULL;
}

/* Return a string corresponding to the error value */
char *xdg_errstr(xdg_error er) {
	switch (er) {
		case xdg_ok:
			return "OK";
		case xdg_alloc:
			return "memory allocation failed";
		case xdg_nohome:
			return "There is no $HOME";
		case xdg_noalluserprofile:
			return "Theres no $ALLUSERSPROFILE is no $ALLUSERSPROFILE";
		case xdg_nopath:
			return "There is no resulting path";
		case xdg_mallformed:
			return "Malfomed path fount";
		default:
			return "unknown";
	}
}


/* ---------------------------------------------------------------- */
#ifdef STANDALONE_TEST
/* test code */

/* Return nz on error */
static int touch(char *name) {
	FILE *fp;

	if ((fp = fopen(name,"w")) == NULL)
		return 1;

	if (fclose(fp))
		return 1;
	 
	return 0;
}

/* Check a file can be opened */
/* Return nz on error */
static int check(char *name) {
	FILE *fp;

	if ((fp = fopen(name,"r")) == NULL)
		return 1;

	if (fclose(fp))
		return 1;
	 
	return 0;
}

/* Delete a path and a file */
/* Return nz on error */
static int delpath(char *path, int depth) {
	int i;
	char *pp;

	for (i = 0; i < depth;) {
//printf("~1 deleting '%s'\n",path);
		if (i == 0) {
			if (unlink(path)) {
//printf("~1 unlink '%s' failed\n",path);
				return 1;
			}
		} else {
			if (rmdir(path)) {
//printf("~1 rmdir '%s' failed\n",path);
				return 1;
			}
		}
		i++;
		if (i == depth)
			return 0;

		if ((pp = strrchr(path, '/')) == NULL)
			return 0;
		*pp = '\000';
	}
	return 0;
}

/* Run a test */
static int runtest(
	xdg_storage_type st,	/* Specify the storage type */
	xdg_scope sc,			/* Scope if write */
	char *pfname,			/* Sub-path and file name */
	char *env,				/* Environment variable being set */
	char *envv,				/* Value to set it to */
	char *defv,				/* default variable needed for read */
	int depth				/* Cleanup depth */
) {
	char *pp;
	xdg_error er;
	char *xval;
	char buf[200];

	if ((xval = getenv(env)) != NULL)		/* Save value before mods */
		xval = strdup(xval);
	if (*env != '\000') {		/* If it is to be set */
		sprintf(buf, "%s=%s",env,envv);
		mputenv(buf);
	}

	printf("\nTesting Variable %s\n",env);
	if ((pp = xdg_bds(&er, st, xdg_write, sc, pfname)) == NULL) {
		printf("Write test failed with %s\n",xdg_errstr(er));
		return 1;
	}
	printf("Create %s %s returned '%s'\n",
		st == xdg_data ? "Data" : st == xdg_data ? "Conf" : "Cache",
		sc == xdg_data ? "User" : "Local", pp);
	if (touch(pp)) {
		printf("Creating file %s failed\n",pp);
		return 1;
	}
	if (check(pp)) {
		printf("Checking file %s failed\n",pp);
		return 1;
	}

	if (sc == xdg_local && *env != '\000') {	/* Add another path */
		sprintf(buf, "%s=xdgtestXXX%c%s",env,SSEP,envv);
		mputenv(buf);
	}

	if (defv != NULL) {
		sprintf(buf, "%s=xdg_NOT_%s",defv,defv);
		mputenv(buf);
	}
	free(pp);

	if ((pp = xdg_bds(&er, st, xdg_read, sc, pfname)) == NULL) {
		printf("Read test failed with %s\n",xdg_errstr(er));
		return 1;
	}
	printf("  Read %s %s returned '%s'\n",
		st == xdg_data ? "Data" : st == xdg_data ? "Conf" : "Cache",
		sc == xdg_data ? "User" : "Local",pp);
	if (check(pp)) {
		printf("Checking file %s failed\n",pp);
		return 1;
	}
	if (delpath(pp, depth)) {
		printf("Warning: Deleting file %s failed\n",pp);
	}
	free(pp);

	/* Restore variables value */
	if (xval == NULL)
		sprintf(buf, "%s=",env);
	else
		sprintf(buf, "%s=%s",env,xval);
	mputenv(buf);

	if (defv != NULL) {
		sprintf(buf, "%s=",defv);
		mputenv(buf);
	}

	return 0;
}

typedef struct {
	xdg_storage_type st;	/* Storage type */
	xdg_scope sc;			/* Scope if write */
	char *defv[2];			/* Default variables needed for user & local tests on read */
	char *envn[10];			/* Environment variable name to set */
} testcase;

int
main() {
	char buf1[200], buf2[200];
	int i, j;

#ifdef NT
	testcase cases[5] = {
		{ xdg_data, xdg_user, {"ALLUSERSPROFILE", NULL},
							{ "XDG_DATA_HOME", "APPDATA", "HOME", "USERPROFILE", NULL } },
		{ xdg_conf, xdg_user, {"ALLUSERSPROFILE", NULL},
							{ "XDG_CONF_HOME", "APPDATA", "HOME", "USERPROFILE", NULL } },
		{ xdg_cache, xdg_user, {NULL, NULL},
							{ "XDG_CACHE_HOME", "APPDATA", "HOME", "USERPROFILE", NULL } },
		{ xdg_data, xdg_local, {NULL, "HOME"},
							{ "XDG_DATA_DIRS", "ALLUSERSPROFILE", NULL } },
		{ xdg_conf, xdg_local, {NULL, "HOME"},
							{ "XDG_CONF_DIRS", "ALLUSERSPROFILE", NULL } }
	};
#else	/* Apple, Unix, Default */
	testcase cases[5] = {
		{ xdg_data, xdg_user, {NULL, NULL},
							{ "XDG_DATA_HOME", "HOME", NULL } },
		{ xdg_conf, xdg_user, {NULL, NULL},
							{ "XDG_CONF_HOME", "HOME", NULL } },
		{ xdg_cache, xdg_user, {NULL, NULL},
							{ "XDG_CACHE_HOME", "HOME", NULL } },
		{ xdg_data, xdg_local, {NULL, "HOME"},
							{ "XDG_DATA_DIRS", "", NULL } },
		{ xdg_conf, xdg_local, {NULL, "HOME"},
							{ "XDG_CONF_DIRS", "", NULL } }
	};
#endif

	/* First clear all the environment variables */
	for (i = 0; i < 5; i++) {
		for (j = 0; ;j++) {
			if (cases[i].envn[j] == NULL)
				break;
			sprintf(buf1, "%s=",cases[i].envn[j]);
			mputenv(buf1);
		}
	}

	/* Then run all the tests */
	for (i = 0; i < 5; i++) {
		for (j = 0; ;j++) {
			if (cases[i].envn[j] == NULL)
				break;
			sprintf(buf1, "xdgtest%d",i);
			sprintf(buf2, "application/%s",cases[i].st == xdg_data ? "data" :
			                   cases[i].st == xdg_conf ? "config" : "cache");
			if (runtest(cases[i].st, cases[i].sc, buf2, cases[i].envn[j],buf1,
			   cases[i].defv[cases[i].sc == xdg_user ? 0 : 1],9))
				exit(1);
		}
	}

	printf("Test completed OK\n");

	return 0;
}

#endif /* STANDALONE_TEST */
