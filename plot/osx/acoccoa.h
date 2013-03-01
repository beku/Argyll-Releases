
#ifndef ACOCCOA_H

/* OS X Coccoa support code for Argyll. */
/* In some places we really prefer not to have OS X code */
/* in separate .m files, so we'll do it all from C. */

#include <Carbon/Carbon.h>

#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 1050
# include <objc/runtime.h>
# include <objc/message.h>
#else
# include <objc/objc-runtime.h>

/* Objective-C runtime compatibility functions for < 10.5 */

/* Create a class definition, but don't register it */
Class CreateClassDefinition(const char *name, const char *superclassName) {
	struct objc_class *meta_class;
	struct objc_class *super_class;
	struct objc_class *new_class;
	struct objc_class *root_class;

	// Ensure that the superclass exists and that someone
	// hasn't already implemented a class with the same name
	//
	super_class = (struct objc_class *)objc_lookUpClass(superclassName);
	if (super_class == nil) {
		printf("failed to lookup '%s'\n",superclassName);
		return Nil;
	}

	if (objc_lookUpClass(name) != nil) {
		return Nil;
	}

	// Find the root class
	root_class = super_class;
	while(root_class->super_class != nil)
		root_class = root_class->super_class;

	// Allocate space for the class and its metaclass
	if ((new_class = calloc(2, sizeof(struct objc_class))) == NULL) {
		return Nil;
	}
	meta_class = &new_class[1];

	// setup class
	new_class->isa = meta_class;
	new_class->info = CLS_CLASS;
	meta_class->info = CLS_META;

	// Create a copy of the class name.
	// For efficiency, we have the metaclass and the class itself
	// to share this copy of the name, but this is not a requirement
	// imposed by the runtime.
	if ((new_class->name = strdup(name)) == NULL) {
		free(new_class);
	}
	meta_class->name = new_class->name;

	// Allocate empty method lists.
	// We can add methods later.
	if ((new_class->methodLists = calloc( 1, sizeof(struct objc_method_list *))) == NULL) {
		free((void *)new_class->name);
		free(new_class);
		return Nil;
	}
	*new_class->methodLists = (struct objc_method_list *) -1;
	if ((meta_class->methodLists = calloc(1, sizeof(struct objc_method_list *))) == NULL) {
		free(new_class->methodLists);
		free((void *)new_class->name);
		free(new_class);
		return Nil;
	}
	*meta_class->methodLists = (struct objc_method_list *) -1;

	// Connect the class definition to the class hierarchy:
	// Connect the class to the superclass.
	// Connect the metaclass to the metaclass of the superclass.
	// Connect the metaclass of the metaclass to the metaclass of the root class.
	new_class->super_class = super_class;
	meta_class->super_class = super_class->isa;
	meta_class->isa = (void *)root_class->isa;

	// Set the sizes of the class and the metaclass.
	new_class->instance_size = super_class->instance_size;
	meta_class->instance_size = meta_class->super_class->instance_size;

	return new_class;
}

/* Add an array of methods. Null terminated by name array */
/* We assume that the class is being created, and that there are */
/* no existing methods. */
BOOL registerDynamicMethods(Class cls, const char *mnames[], IMP mimps[], const char *mtypes[]) {
	int i, nmeth;
	struct objc_method_list *methodList;
   
	/* Count the number of methods */
	for (nmeth = 0; mnames[nmeth] != NULL && mnames[nmeth][0] != '\000'; nmeth++)
		;

	/* Allocate an array */
	methodList = malloc(sizeof(struct objc_method_list) + (nmeth-1) * sizeof(struct objc_method));
       
	methodList->method_count = nmeth;
	for (i = 0; i < nmeth; i++) {
		// Get or register the selector for the method name
		SEL methodSEL = SELUID(mnames[i]);
	
		// Registering the method seems to register the selector
		if (ISSELECTOR(methodSEL) == NO) {
			methodSEL = sel_registerName(mnames[i]);
		}
  	 
		// Fill out the method list
		methodList->method_list[i].method_name = methodSEL;
		methodList->method_list[i].method_imp = mimps[i];
		methodList->method_list[i].method_types = strdup(mtypes[i]);
	}
	
	// Register our methods
	class_addMethods((Class)cls, methodList);

	return YES;
}

/* Add an array of instance variables. Null terminated by name array */
/* We assume that the class is being created, and that there are */
/* no existing methods. */
BOOL registerDynamicVariables(Class cls, const char *names[], size_t sizes[], const char *types[]) {
	int i, nvar = 1;
	int vsize;
	struct objc_ivar *ivarp;

	/* Count the number of variables */
	for (nvar = 0; names[nvar] != NULL && names[nvar][0] != '\000'; nvar++)
		;

	vsize = sizeof(struct objc_ivar_list) + (nvar - 1) * sizeof(struct objc_ivar);
	cls->ivars = calloc(vsize, 1);
	cls->ivars->ivar_count = nvar;

	for (i = 0; i < nvar; i++) {
		int abytes;
		ivarp = &cls->ivars->ivar_list[i];
	
		/* Set the variable information */
		ivarp->ivar_name = strdup(names[i]);
		ivarp->ivar_type = strdup(types[i]);
	
		/* Align the offset for this variable to it's size, limiting to 64 bits */
		if ((abytes = sizes[i]) > 8)
			abytes = 8;
		cls->instance_size = (cls->instance_size + abytes-1) & ~(abytes-1);
		ivarp->ivar_offset = (int)cls->instance_size;
		cls->instance_size += sizes[i];
	}

    return YES;
}

#endif /* __MAC_OS_X_VERSION_MAX_ALLOWED >= 1050 */

extern id NSApp;

/* Create a class */
BOOL registerClass(
const char *name,			/* Name of class being created */
const char *supername,		/* Name of superclass */
const char *methnames[],	/* List of method names, empty string terminated */
IMP methimps[],				/* Method implementations */
const char *methsigs[],			/* Method signatures */
const char *varnames[],		/* List of variable names, empty string terminated */
size_t varsizes[],			/* Variable size in bytes */
const char *varsigs[]		/* Variable signatures */
) {
	int i;
	Class nclass;

#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 1050
	nclass = objc_allocateClassPair((Class)objc_getClass(supername), name, 0);
	if (nclass == Nil) {
		return NO;
	}

	for (i = 0; methnames[i] != NULL && methnames[i][0] != '\000'; i++) {
		class_addMethod(nclass, sel_getUid(methnames[i]), methimps[i], methsigs[i]);
	}

	// Should check return value is YES
	for (i = 0; varnames[i] != NULL && varnames[i][0] != '\000'; i++) {
		int asize;
		while( (1 << asize) < varsizes[i] && asize < 8)
			asize++;
		class_addIvar(nclass, varnames[i],  varsizes[i], asize, varsigs[i]); 
	}

	// Must be called after adding instance variable
	objc_registerClassPair(nclass);
#else
	/* Use compat.h functions to do the dirty work */
	// Should check the return value!
	if ((nclass = CreateClassDefinition(name, supername)) == Nil) {
		return NO;
	}

	registerDynamicMethods(nclass, methnames, methimps, methsigs);

	registerDynamicVariables(nclass, varnames, varsizes, varsigs);

	// Register the class with the runtime.
	objc_addClass(nclass);
#endif
	return YES;
}

/* ------------------------------------------------ */
/* One of the primary disadvantages of coding Coccoa in C */
/* is the lack of compatible .h files. We have to make our own.. */

/* ------------------------------------------------ */

/* Foundation stuff */

#ifndef __OBJC__

#if !defined(FSTATIC_INLINE)
# if defined (__GNUC__) && (__GNUC__ == 4)
#  define FSTATIC_INLINE static __inline__ __attribute__((always_inline))
# else
#  define FSTATIC_INLINE static __inline__
# endif
#endif

#ifdef __LP64__
typedef double NSFloat;
#else
typedef float NSFloat;
#endif

typedef struct _NSPoint {
    NSFloat x;
    NSFloat y;
} NSPoint;

FSTATIC_INLINE NSPoint NSMakePoint(NSFloat x, NSFloat y) {
    NSPoint r;
    r.x = x;
    r.y = y;
    return r;
}

typedef struct _NSSize {
    NSFloat width;
    NSFloat height;
} NSSize;

FSTATIC_INLINE NSSize NSMakeSize(NSFloat w, NSFloat h) {
    NSSize r;
    r.width = w;
    r.height = h;
    return r;
}


typedef struct _NSRect {
    NSPoint origin;
    NSSize size;
} NSRect;

FSTATIC_INLINE NSRect NSMakeRect(NSFloat x, NSFloat y, NSFloat w, NSFloat h) {
    NSRect r;
    r.origin.x = x;
    r.origin.y = y;
    r.size.width = w;
    r.size.height = h;
    return r;
}

#endif /* !__OBJC__ */

/* ------------------------------------------------ */
/* Constats for NSString class */

/* 10.4 and latter */
typedef enum {
	NSStringDrawingUsesLineFragmentOrigin        = (1 << 0),
	NSStringDrawingUsesFontLeading               = (1 << 1),
	NSStringDrawingDisableScreenFontSubstitution = (1 << 2),
	NSStringDrawingUsesDeviceMetrics             = (1 << 3),
	NSStringDrawingOneShot                       = (1 << 4),
	NSStringDrawingTruncatesLastVisibleLine      = (1 << 5)
} NSStringDrawingOptions;

/* ------------------------------------------------ */
/* Constats for NSApplication class */

typedef enum {
	NSApplicationPresentationDefault                    = 0,
	NSApplicationPresentationAutoHideDock               = (1 <<  0),
	NSApplicationPresentationHideDock                   = (1 <<  1),
	NSApplicationPresentationAutoHideMenuBar            = (1 <<  2),
	NSApplicationPresentationHideMenuBar                = (1 <<  3),
	NSApplicationPresentationDisableAppleMenu           = (1 <<  4),
	NSApplicationPresentationDisableProcessSwitching    = (1 <<  5),
	NSApplicationPresentationDisableForceQuit           = (1 <<  6),
	NSApplicationPresentationDisableSessionTermination  = (1 <<  7),
	NSApplicationPresentationDisableHideApplication     = (1 <<  8),
	NSApplicationPresentationDisableMenuBarTransparency = (1 <<  9),
	NSApplicationPresentationFullScreen                 = (1 << 10),
	NSApplicationPresentationAutoHideToolbar            = (1 << 11)
} NSApplicationPresentationOptions;

typedef enum {
	NSTerminateCancel = 0,
	NSTerminateNow = 1,
	NSTerminateLater = 2
} NSApplicationTerminateReply;

/* ------------------------------------------------ */
/* Constats for NSWindow class */

enum {
    NSBorderlessWindowMask               = 0,
    NSTitledWindowMask                   = 1 << 0,
    NSClosableWindowMask                 = 1 << 1,
    NSMiniaturizableWindowMask           = 1 << 2,
    NSResizableWindowMask                = 1 << 3,
	NSTexturedBackgroundWindowMask       = 1 << 8
#if MAC_OS_X_VERSION_MAX_ALLOWED >= MAC_OS_X_VERSION_10_4
	,NSUnscaledWindowMask                = 1 << 11,
	 NSUnifiedTitleAndToolbarWindowMask  = 1 << 12
#endif
};

/* types of window backing store */
typedef enum {
    NSBackingStoreRetained       = 0,
    NSBackingStoreNonretained    = 1,
    NSBackingStoreBuffered       = 2
} NSBackingStoreType;


/* ------------------------------------------------ */
/* Convenience functions */

/* Transform process to interact with descktop */
void transProcess() {
	OSStatus stat;
	ProcessSerialNumber psn = { 0, 0 };

	if (GetCurrentProcess(&psn) != noErr) {
		/* Transform the process so that the desktop interacts with it properly. */
		/* We don't need resources or a bundle if we do this. */
		if (psn.lowLongOfPSN != 0 && (stat = TransformProcessType(&psn,
			               kProcessTransformToForegroundApplication)) != noErr)
			fprintf(stderr,"TransformProcess failed with code %d\n",stat);
	}
}

/* Send a message to the object */
#define sendMsg(oid, msg, ...) objc_msgSend(oid, sel_getUid(msg), ##__VA_ARGS__)

/* alloc and init a new object */
id newObject(const char *cname) {
	id rv;
	rv = objc_msgSend(objc_getClass(cname), sel_getUid("alloc"));
	rv = objc_msgSend(rv, sel_getUid("init"));
	return rv;
}
/* release an object */
void delObject(id oid) {
	objc_msgSend(oid, sel_getUid("release"));
}

/* dealloc super */
void delObjectSuper(id oid) {
	struct objc_super ss;
	ss.receiver = oid;
#ifdef __OBJC__
	ss.super_class = sendMsg(oid, "superclass");
#else
	ss.class = (Class)sendMsg(oid, "superclass");
#endif
	objc_msgSendSuper(&ss, sel_getUid("dealloc"));
}

/* Create an NSString from a C string */
id newNSString(const char *cstr) {
	id str;

	str = objc_msgSend(objc_getClass("NSString"), sel_getUid("alloc"));
	str = objc_msgSend(str, sel_getUid("initWithUTF8String:"), cstr);

	return str;
}

#define ACOCCOA_H
#endif /* ACOCCOA_H */
