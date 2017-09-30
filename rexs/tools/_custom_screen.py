import curses
import curses.textpad
import curses.ascii
import curses.panel
from collections import defaultdict
def validator(ch):
    if ch==curses.ascii.NL:
        return curses.ascii.BEL
    if not curses.ascii.isprint(ch):
        reprint()
    return ch


class Screen:
    def __init__(self, welcome=""):
        self.S = curses.initscr()
        curses.cbreak()
        curses.noecho()
        self.S.keypad(1)
        self.offsety = 2
        self._y, self._x = self.S.getmaxyx()
        if welcome:
            self.printlines(welcome)
            self.S.addstr(self._y-2, 2, "Hit any key...", curses.A_BOLD)
            self.S.border(0)
            self.S.getch()
    def __enter__(self):
        return self
    def __exit__(self,a,b,c):
        curses.nocbreak()
        self.S.keypad(0)
        curses.echo()
        curses.endwin()
    def printlines(self, lines, offsetx=2, offsety=None, fmt=0):
        if isinstance(lines, str):
            lines = [lines]
        if offsety is None:
            offsety = self.offsety
        line = 0
        width = self._x - offsetx
        for i, msg in enumerate(lines):
            #self.S.addnstr(offsety + line, offsetx - 1, "-", 1)
            if not ((offsety + line) < self._y) or not ((offsetx+1) < self._x):
                break
            while msg:
                try:
                    self.S.addnstr(offsety + line, offsetx, msg, width, fmt)
                except:
                    break
                msg = msg[width:]
                line+=1
        self.offsety = offsety + line
        self.S.refresh()
        return len(lines) - (i+1)
    def prompt(self, question, y=2, x=2, dy=3, dx=None, options=()):
        if dx is None:
            dx = self._x - x - 2
        def reprint():
            maxlen=0
            for option in options:
                self.printlines(option, x+maxlen+4, y+dy+3)
                maxlen += max(map(len, option))
            self.S.addstr(y, x, question)
            curses.textpad.rectangle(self.S, y+1, x, y+dy+2, x+dx)
        reprint()
        self._reprint = reprint
        self.S.refresh()
        self.win = curses.newwin(dy, dx-1, y+2, x+1)
        self.tb = curses.textpad.Textbox(self.win)
        text = self.tb.edit(validator)
        del self.tb
        del self._reprint
        #inp = self.S.getstr(y+dy, x+dx)
        return text
    
    def menu(self, items, position=0, info={}, showall=False, selected=None,
                   readonly = [], submit={}, **kwargs):
        """
            Opens a list-like menu:
            
            Inputs:
                items : list
                    list of items to choose from
                position : int
                    starting position of the curser
                info : dict
                    dictionary containing information for a set of items
                selected : list
                    list of marked items, will be modified during call
                readonly : list
                    list of items that cannot be selected
                submit : list or dict
                    list or dict of extra items that contain the submit options
                    if dict, the values contain the keys that cause the 
                    corresponding submission when pressed
                
        """
        if isinstance(submit, list):
            submit = dict.fromkeys(submit)
        if position in items:
            position = items.index(position)
        elif not isinstance(position, int):
            position = 0
        submit = defaultdict(str, submit)
        submit["cancel"] = 27
        mark = 8*"="
        items.append(mark)
        readonly.append(mark)
        items+=submit.keys()
        buttons = dict([(v,k) for (k,v) in submit.iteritems()])
        key_unselect = kwargs.get("unselect", ord("-"))
        key_groupselect = kwargs.get("groupselect", ord("+"))
        key_groupunselect = kwargs.get("groupunselect", ord("_"))
        #info["cancel"] = str(submit)
        self.win = self.S.subwin(0,0)
        self.win.keypad(1)
        self.panel = curses.panel.new_panel(self.win)
        #self.panel.hide()
        curses.panel.update_panels()
        self.panel.top()
        self.panel.show()
        self.win.clear()
        values = list(items) # make a copy
        items = map(str, items)
        maxlen = max(map(len, items)) + 1
        space1 = len("%d. "%len(items)) + 1
        space2 = kwargs.get("space2", 3)
        infolen = self._x - 1 - maxlen - space2
        def navigate(n, pos):
            pos += n
            if pos < 0:
                pos = 0
            elif pos >= len(items):
                pos = len(items)-1
            return pos
        
        digits = [ord(str(i)) for i in xrange(min(len(items), 9))]
        
        while True:
            self.win.clear()
            self.win.refresh()
            curses.doupdate()
            self.draw_title()
            #if selected!=None:
            #    self.win.addstr(0, 0, 
            #        "Use: <space> - mark item/variable | <enter> - select |"\
            #        " q - quit",
            #        curses.A_UNDERLINE)
            #start = max(0, position - (self._y-7))
            start = max(0, position - (self._y-self._y/2))
            start = min(start, len(items) - (self._y - 5))
            #print start
            for index, item in enumerate(items):
                if index == position:
                    mode = curses.A_REVERSE
                else:
                    mode = curses.A_NORMAL 
                if index+1-start>=self._y-4:
                    break
                if index<start:
                    continue
                if item not in (list(submit) + [mark]):
                    self.win.addstr(1+index-start, 1, '%d.'%index, mode)
                self.win.addstr(1+index-start, space1, '%s'%item, mode)
                if selected!=None and item in selected:
                    self.win.addstr(1+index-start, 0, "*", curses.A_BOLD)
                if item in info and showall:
                    iteminfo = info[item]
                    self.win.addstr(1+index-start, space1 + space2 + maxlen, 
                                    iteminfo[:infolen], mode)
            if items[position] in info:
                iteminfo = info[items[position]]
                self.win.addstr(self._y-3, 0, iteminfo, curses.A_BOLD)

            key = self.win.getch()
            if key == curses.KEY_UP:
                position = navigate(-1, position)

            elif key == curses.KEY_DOWN:
                position = navigate( 1, position)
            elif key == curses.KEY_PPAGE:
                position = navigate( -self._y//2, position)
            elif key == curses.KEY_NPAGE:
                position = navigate( self._y//2, position)
            elif key == curses.KEY_END:
                position = len(items)-1
            elif key == curses.KEY_HOME:
                position = 0
            elif key in digits:
                position = digits.index(key)
            elif key == key_unselect:
                del selected[:]
            elif key == key_groupselect:
                for p in xrange(position, 0, -1):
                    if values[p] in selected:
                        break
                    else:
                        selected.append(values[p])
            elif key == key_groupunselect:
                for p in xrange(position, 0, -1):
                    if values[p] in selected:
                        selected.remove(values[p])
                    else:
                        break
            elif key in [curses.KEY_ENTER, ord('\n')]:
                if selected and values[position] not in submit:
                    return "__defaultaction"
                else:
                    break
            elif key in submit.values():
                action = buttons[key]
                return action
            elif values[position] in readonly:
                continue
            elif key in [ord(' ')] and selected!=None:
                if values[position] in submit:
                    continue
                elif values[position] in selected:
                    selected.remove(values[position])
                else:
                    selected.append(values[position])
        
        self.win.clear()
        self.panel.hide()
        curses.panel.update_panels()
        curses.doupdate()
        return values[position]
    
    def draw_title(self, title=None):
        if title is None and hasattr(self, "title"):
            title = self.title
        elif title is None:
            return
        else:
            self.title = title
        self.S.addstr(0,0, title, curses.A_BOLD + curses.A_UNDERLINE)
    
    def resize(self, signum=None, frame=None):
        self.clear()
        curses.endwin()
        self.S = curses.initscr()
        self._y, self._x = self.S.getmaxyx()
        #curses.resizeterm(self._y, self._x)
        #self.S.resize(self._y, self._x)
        #self.S.addstr(9,5, str(resize))
        winsize = "rows: %i, cols: %i"%(self._y, self._x)
        self.S.addstr(self._y-1, self._x-len(winsize)-1, winsize)
        #self.S.addstr(0, self._x-len(winsize)-1, winsize)
        if hasattr(self, "_reprint"):
            self._reprint()
        if hasattr(self, "tb"):
            text = self.tb.gather()
            [self.tb.do_command(curses.ascii.BS) for ch in text]
            text = text.strip()
            [self.tb.do_command(ch) for ch in text]
        self.draw_title()
        self.S.refresh()
    def clear(self):
        self.S.clear()
        self.offsety = 2


