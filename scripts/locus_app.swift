import Cocoa

// Locus.app — a minimal native Cocoa app: shows in the Dock, never bounces, has a real
// Quit (Cmd-Q / menu / Dock), starts the Locus server on launch, opens the browser, and
// stops the server when you quit.
//
// @REPO@ and @UV@ are substituted by scripts/build_macos_app.sh at build time.

let REPO = "@REPO@"
let UV = "@UV@"

func stopServer() {
    let pk = Process()
    pk.executableURL = URL(fileURLWithPath: "/usr/bin/pkill")
    pk.arguments = ["-f", "locus serve api"]
    try? pk.run()
    pk.waitUntilExit()
}

final class AppDelegate: NSObject, NSApplicationDelegate {
    let task = Process()
    let url = URL(string: "http://127.0.0.1:8787")!

    func applicationDidFinishLaunching(_ notification: Notification) {
        startServer()
        DispatchQueue.main.asyncAfter(deadline: .now() + 2.5) { [weak self] in
            if let u = self?.url { NSWorkspace.shared.open(u) }
        }
    }

    func startServer() {
        let reports = REPO + "/data/reports"
        try? FileManager.default.createDirectory(atPath: reports, withIntermediateDirectories: true)
        let logPath = reports + "/serve.log"
        if !FileManager.default.fileExists(atPath: logPath) {
            FileManager.default.createFile(atPath: logPath, contents: nil)
        }
        task.executableURL = URL(fileURLWithPath: UV)
        task.arguments = ["run", "locus", "serve", "api"]
        task.currentDirectoryURL = URL(fileURLWithPath: REPO)
        if let fh = FileHandle(forWritingAtPath: logPath) {
            fh.seekToEndOfFile()
            task.standardOutput = fh
            task.standardError = fh
        }
        do { try task.run() } catch { NSLog("Locus: failed to start server: \(error.localizedDescription)") }
    }

    // Clicking the Dock icon while running re-opens the page.
    func applicationShouldHandleReopen(_ sender: NSApplication, hasVisibleWindows flag: Bool) -> Bool {
        NSWorkspace.shared.open(url)
        return true
    }

    func applicationWillTerminate(_ notification: Notification) {
        stopServer()
    }
}

let app = NSApplication.shared
let delegate = AppDelegate()
app.delegate = delegate
app.setActivationPolicy(.regular)   // Dock icon + app menu

// Minimal menu so Cmd-Q / "Quit Locus" works.
let mainMenu = NSMenu()
let appItem = NSMenuItem()
mainMenu.addItem(appItem)
let appMenu = NSMenu()
appMenu.addItem(withTitle: "Quit Locus", action: #selector(NSApplication.terminate(_:)), keyEquivalent: "q")
appItem.submenu = appMenu
app.mainMenu = mainMenu

// Also stop the server on SIGTERM/SIGINT (covers kill/logout). Force-quit (SIGKILL) can't be caught.
signal(SIGTERM, SIG_IGN)
signal(SIGINT, SIG_IGN)
var sigSources: [DispatchSourceSignal] = []
for s in [SIGTERM, SIGINT] {
    let src = DispatchSource.makeSignalSource(signal: s, queue: .main)
    src.setEventHandler { stopServer(); exit(0) }
    src.resume()
    sigSources.append(src)
}

app.run()
