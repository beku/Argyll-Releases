
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __OBJC__
# include <Foundation/Foundation.h>
# include <AppKit/AppKit.h>
#endif

#include "acocoa.h"		/* Argyll Cocoa support functions and defines */

/*
- (void) someMethod:...
{
    va_list va;

    va_start(va, _cmd);

    // process all args with va_arg

    va_end(va);
}
 */

/* Our static instance variables for AppDelegate */
typedef struct {
    id window;		/* NSWindow */
    id view;		/* NSView */
} cntx_t;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void MyView_setCntx(id self, SEL _cmd, void *val) {
	object_setInstanceVariable(self, "cntx", val);
}

void MyView_drawRect(id self, SEL _cmd, NSRect rect) {
	double w = rect.size.width, h = rect.size.height;

	id aPath = sendClassMsg("NSBezierPath", "bezierPath");
	sendMsg(aPath, "setLineWidth:", 2.0);
	sendMsg(aPath, "moveToPoint:", NSMakePoint(0.0, 0.0));
	sendMsg(aPath, "lineToPoint:", NSMakePoint(0.9 * w, 0.9 * h));
	sendMsg(aPath, "appendBezierPathWithRect:", NSMakeRect(0.5 * w, 0.5 * h,
		                                                   0.7 * w, 0.6 * h));
	sendMsg(aPath, "stroke");

	{
		id att = newObject("NSDictionary");
		id str = newNSString("String");
		sendMsg(str, "drawAtPoint:withAttributes:", NSMakePoint(0.1 * w, 0.1 * h), att); 
	}
}

/* Clean up */
void MyView_dealloc(id self, SEL _cmd) {
	delObjectSuper(self);
}

// Create our custom NSView */
void createMyView() {
	method_info minfo[] = {
		{ "setCntx:",          (IMP)MyView_setCntx,          "v@:^v" },
		{ "drawRect:",         (IMP)MyView_drawRect,         "v@:@"  },
		{ "dealloc",           (IMP)MyView_dealloc,          "v@:"   },
		{ "" }
	};

	variable_info vinfo[] = {
		{ "cntx", sizeof(void *), "^v" },
		{ "" }
	};

	registerClass("MyView", "NSView", minfo, vinfo);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void MyWin_setCntx(id self, SEL _cmd, void *val) {
	object_setInstanceVariable(self, "cntx", val);
}

void MyWin_keyDown(id self, SEL _cmd, id event) {
	int etype;
	id nresp;
	id str;
	const char *cstr;

	etype = (int)sendMsg(event, "type");
	str = sendMsg(event, "characters");
	cstr = cNSString(str); 
	printf("Got Window KeyDown type %d, chars '%s'\n",etype, cstr);

	if (cstr[0] == ' ')
		sendMsg(NSApp, "terminate:", self);
}

/* Clean up */
void MyWin_dealloc(id self, SEL _cmd) {
	delObjectSuper(self);
}

// Create our custom NSWin */
void createMyWin() {
	method_info minfo[] = {
		{ "setCntx:",             (IMP)MyWin_setCntx,              "v@:^v", },
		{ "keyDown:",             (IMP)MyWin_keyDown,              "v@:@",  },
		{ "dealloc",              (IMP)MyWin_dealloc,              "v@:"    },
		{ "" }
	};

	variable_info vinfo[] = {
		{ "cntx", sizeof(void *), "^v" },
		{ "" }
	};

	registerClass("MyWin", "NSWindow", minfo, vinfo);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Create all the bits */
void AppDel_willFinishLaunching(id self, SEL _cmd, id notification) {
	cntx_t *cx;
	id label;		/* NSTextField */

	cx = calloc(1, sizeof(cntx_t));

	// Set cntx to to allocated structure
	object_setInstanceVariable(self, "cntx", (void *)cx);

	/* Create Window */
    cx->window = sendClassMsg("MyWin", "alloc");
    cx->window = sendMsg(cx->window,
	                             "initWithContentRect:styleMask:backing:defer:",
	                             NSMakeRect(300, 300, 200, 100),
                                 NSTitledWindowMask
	                           | NSClosableWindowMask
	                           | NSMiniaturizableWindowMask
	                           | NSResizableWindowMask,
                                 NSBackingStoreBuffered,
	                             YES);

	/* Make the background white */
	sendMsg(cx->window, "setBackgroundColor:", sendClassMsg("NSColor","whiteColor"));

	/* Add title */
	sendMsg(cx->window, "setTitle:", newNSString("Hello World"));
	
#ifdef NEVER
	/* Create Label */
    label = sendClassMsg("NSTextField", "alloc");
    label = sendMsg(label, "initWithFrame:", NSMakeRect(30, 30, 80, 30));

    sendMsg(label, "setSelectable:", NO);
    sendMsg(label, "setBezeled:", NO);
    sendMsg(label, "setDrawsBackground:", NO);
    sendMsg(label, "setStringValue:", newNSString("Hello World"));

	/* Hmm. How does this work ? */
    cx->view = sendMsg(cx->window, "contentView");
    sendMsg(cx->view, "addSubview:", label);
#else
	/* Use our custom view to draw contents */
    cx->view = newObject("MyView");

    sendMsg(cx->view, "setCntx:", (void *)cx);
    sendMsg(cx->window, "setContentView:", cx->view);

//	sendMsg(cx->window, "setInitialFirstResponder:", cx->view);
#endif

	// Window methods:

	// sendEvent: gets messages.

	// Set above the screen saver
	// [aWindow setLevel:NSScreenSaverWindowLevel + 1];

	// setCollectionBehavior: NSWindowCollectionBehaviorIgnoresCycle
	//                        NSWindowCollectionBehaviorFullScreenPrimary
	//                        NSWindowCollectionBehaviorStationary				
	//                        NSWindowCollectionBehaviorIgnoresCycle				

	// center
	// mouseDownCanMoveWindow 
	// setMovable: 
	// setContentMinSize: and setContentMaxSize:
	//
    // NSRect frame = [myWindow frame];
    // if (frame.size.width <= MIN_WIDTH_WITH_ADDITIONS)
    //     frame.size.width = MIN_WIDTH_WITH_ADDITIONS;
    // frame.size.height += ADDITIONS_HEIGHT;
    // frame.origin.y -= ADDITIONS_HEIGHT;
    // [myWindow setFrame:frame display:YES animate:YES];
    // objc_msgSend(label, "release"));
	//
	// setExcludedFromWindowsMenu:YES
	//
	// setBackgroundColor: and setAlphaValue:

	// WindowImage object:
	// setDepthLimit:

	// Setting cursor:
	// NSTrackingArea class, along with the cursorUpdate: method of the NSResponder class
	
}

/* Map the window */
void AppDel_didFinishLaunching(id self, SEL _cmd, id notification) {
    cntx_t *cx;

	object_getInstanceVariable(self, "cntx", (void **)&cx);
    sendMsg(cx->window, "makeKeyAndOrderFront:", self);

#ifdef NEVER		/* Test terminate */
	sleep(5);
	sendMsg(NSApp, "terminate:", self);
#endif
}

/* Should the application terminate ? */
NSApplicationTerminateReply AppDel_shouldTerminate(id self, SEL _cmd, id notification) {
	return NSTerminateNow;
}

/* Clean up */
void AppDel_dealloc(id self, SEL _cmd) {
    cntx_t *cx;

	object_getInstanceVariable(self, "cntx", (void **)&cx);
    delObject(cx->window);
	delObjectSuper(self);
}

// Create the delegate class that implements our window
void createAppDelClass() {
	method_info minfo[] = {
		{ "applicationWillFinishLaunching:", (IMP)AppDel_willFinishLaunching, "v@:@" },
		{ "applicationDidFinishLaunching:",  (IMP)AppDel_didFinishLaunching,  "v@:@" },
		{ "applicationShouldTerminate:",     (IMP)AppDel_shouldTerminate,     "i@:@" },
		{ "dealloc",                         (IMP)AppDel_dealloc,             "v@:"  },
		{ "" }
	};

	variable_info vinfo[] = {
		{ "cntx", sizeof(void *), "^v" },
		{ "" }
	};

	registerClass("AppDelegate", "NSObject", minfo, vinfo);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - */
int main(int argc, char** argv) {
	id pool;
	id appDelObj;

	/* Transform process so that it interacts with desktop properly */
	transProcess();

	/* Create an autorelease pool */
	pool = newObject("NSAutoreleasePool");

	// Create all the classes we override
	createAppDelClass();
	createMyWin();
	createMyView();

	// Get our shared NSApplication and start it
	sendClassMsg("NSApplication", "sharedApplication");

	if (NSApp == NULL) {
		fprintf(stderr,"Failed to initialized NSApplication...  terminating...\n");
		return -1;
	}

	/* Set a delegate to create the window */
	appDelObj = newObject("AppDelegate");
	sendMsg(NSApp, "setDelegate:", appDelObj);

// if running on 10.7:
//	sendMsg(NSApp, "disableRelaunchOnLogin"));

	/* Call the run loop */
	sendMsg(NSApp, "run");

	// detachDrawingThread:toTarget:withObject:

	/* To terminate:
	sendMsg(NSApp, "terminate:", self);
	*/

	/* We're done with the pool */
    delObject(pool);

	return EXIT_SUCCESS;
}

#ifdef NEVER

- (void)drawRect:(NSRect)rect
{
  NSRect r = NSMakeRect(10, 10, 50, 60);
  NSBezierPath *bp = [NSBezierPath bezierPathWithRect:r];
  NSColor *color = [NSColor blueColor];
  [color set];
  [bp fill];
}

#endif
