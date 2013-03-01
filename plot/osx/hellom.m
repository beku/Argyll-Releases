
#include <Foundation/Foundation.h>
#include <AppKit/AppKit.h>

#include <objc/objc-runtime.h>		/* For diagnostics */

typedef struct {
	id window;		/* NSWindow */
	id view;		/* NSTextField */
} cntx_t;

// - - - - - - - - - - - - - - - - - - - - - - - - - 
@interface MyView : NSView {
	void *cntx;
}
- (void)setCntx:(void *)cntx;
@end

@implementation MyView

- (void)setCntx:(void *)val {
	cntx = val;
}

- (void)drawRect:(NSRect)rect {
//	NSGraphicsContext* aContext = [NSGraphicsContext currentContext];

	[[NSColor redColor] setStroke];

	NSBezierPath* aPath = [NSBezierPath bezierPath];
	[aPath setLineWidth:2.0];
	[aPath moveToPoint:NSMakePoint(0.0, 0.0)];
	[aPath lineToPoint:NSMakePoint(100.0, 100.0)];
	[aPath appendBezierPathWithRect:NSMakeRect(20.0, 160.0, 80.0, 50.0)];

	[aPath stroke];

	NSDictionary *att = [NSDictionary new];

	[ @"String" drawAtPoint: NSMakePoint(10.0, 10.0) withAttributes: att ];

/* 

[@"Hello" drawInRect:r withAttributes:[NSDictionary
dictionaryWithObjectsAndKeys:
[NSColor redColor], NSForegroundColorAttributeName,
[NSFont systemFontOfSize:24], NSFontAttributeName,
nil]];

*/
}

@end

/* To trigger an update:

Send a setNeedsDisplayInRect: or setNeedsDisplay: message to the
view. Sending either of these messages marks part or all of the view as invalid

 */

// - - - - - - - - - - - - - - - - - - - - - - - - - 

@interface MyWin : NSWindow {
	void *cntx;
}
- (void)setCntx:(void *)cntx;
@end

@implementation MyWin

- (void)setCntx:(void *)val {
	cntx = val;
}

- (void)keyDown:(NSEvent *)event {
	printf("Got Window KeyDown type %d char %s\n",[event type],
	    [[event characters] cStringUsingEncoding:NSASCIIStringEncoding]);
}

- (BOOL)windowShouldClose:(id)sender {
	printf("Got Window windowShouldClose\n"); 
	[NSApp terminate: nil];
	return YES;
}

@end

// - - - - - - - - - - - - - - - - - - - - - - - - - 

@interface AppDelegate : NSObject {
	void *cntx;
}
@end

@implementation AppDelegate
- (void) applicationWillFinishLaunching: (NSNotification *)not {
	cntx_t *cx;

	cx = calloc(1, sizeof(cntx_t));

	cntx = (void *)cx;

	/* Create Window */
	cx->window = [[MyWin alloc] initWithContentRect: NSMakeRect(300, 300, 200, 100)
                                        styleMask: (NSTitledWindowMask |
	                                                NSClosableWindowMask |
                                                    NSMiniaturizableWindowMask |
                                                    NSResizableWindowMask)
                                          backing: NSBackingStoreBuffered
                                            defer: YES];
	[cx->window setTitle: @"Hello World"];

#ifdef NEVER
	/* Create Label */
	cx->label = [[NSTextField alloc] initWithFrame: NSMakeRect(30, 30, 80, 30)]; 
	[cx->label setSelectable: NO];
	[cx->label setBezeled: NO];
	[cx->label setDrawsBackground: NO];
	[cx->label setStringValue: @"Hello World"];

	[[cx->window contentView] addSubview: cx->label];

	[cx->label release];

#else
	cx->view = [MyView new];
	[cx->view setCntx:(void *)cx];
	[cx->window setContentView: cx->view];
#endif
	// [window setContentView:customView]
}

- (void) applicationDidFinishLaunching: (NSNotification *) not {
	cntx_t *cx = (cntx_t *)cntx;
	[cx->window makeKeyAndOrderFront: self];
}

- (void) dealloc {
	cntx_t *cx = (cntx_t *)cntx;
	[cx->window release];
	[super dealloc];
}
@end

// - - - - - - - - - - - - - - - - - - - - - - - - - 

int main (int argc, const char **argv)
{ 
	NSAutoreleasePool *pool;
	id appDelObj;
   
	if (NSApp == nil) {
		OSStatus stat;
		ProcessSerialNumber psn = { 0, 0 };

		if (GetCurrentProcess(&psn) == noErr) {
			/* Transform the process so that the desktop interacts with it properly. */
			/* We don't need resources or a bundle if we do this. */
			if (psn.lowLongOfPSN != 0 && (stat = TransformProcessType(&psn,
				               kProcessTransformToForegroundApplication)) != noErr)
				fprintf(stderr,"TransformProcess failed with code %d\n",stat);
		}

		pool = [NSAutoreleasePool new];

		[NSApplication sharedApplication];
		appDelObj = [AppDelegate new];
		[NSApp setDelegate: appDelObj];

		// Run the event loop until done
		[NSApp run];

		[pool release];
	}

	// To terminate:
	// [NSApp terminate: nil];
}


