"""Background job queue: one worker thread, so only one assembly runs at a time."""
import logging
import queue
import threading
from pathlib import Path

from db import mongo
from pipeline import runner
from pipeline.errors import MESSAGES, TOOL_FAILED

log = logging.getLogger("amr.jobs")


class JobManager:
    def __init__(self):
        self._q: queue.Queue = queue.Queue()
        self._thread: threading.Thread | None = None
        self._lock = threading.Lock()

    def start(self) -> None:
        with self._lock:
            if self._thread and self._thread.is_alive():
                return
            self._thread = threading.Thread(target=self._loop, name="amr-worker", daemon=True)
            self._thread.start()

    def submit(self, isolate_id: str, r1: Path, r2: Path, species: str) -> None:
        self.start()
        self._q.put((isolate_id, r1, r2, species))

    def pending(self) -> int:
        return self._q.qsize()

    def _loop(self) -> None:
        while True:
            isolate_id, r1, r2, species = self._q.get()
            try:
                runner.run_isolate(isolate_id, r1, r2, species)
            except Exception:
                log.exception("worker error on %s", isolate_id)
            finally:
                self._q.task_done()

    def recover(self) -> None:
        """After a restart: re-queue isolates still queued (if their reads exist), fail the interrupted ones."""
        try:
            for iso in mongo.find_isolates_by_status("running"):
                mongo.set_isolate_status(iso["_id"], "failed", error_code=TOOL_FAILED,
                                         error_message=MESSAGES[TOOL_FAILED])
            for iso in mongo.find_isolates_by_status("queued"):
                d = runner.upload_dir(iso["_id"])
                files = sorted(d.glob("R[12].*")) if d.exists() else []
                if len(files) == 2:
                    self.submit(iso["_id"], files[0], files[1], iso.get("species", ""))
                else:
                    mongo.set_isolate_status(iso["_id"], "failed", error_code=TOOL_FAILED,
                                             error_message=MESSAGES[TOOL_FAILED])
        except Exception:
            log.warning("job recovery skipped", exc_info=True)


manager = JobManager()
