import multiprocessing
import sys


def main() -> int:
    multiprocessing.freeze_support()
    if len(sys.argv) > 1:
        from .cli import main as cli_main
        return cli_main(sys.argv[1:])
    from .gui import main as gui_main
    return gui_main()


if __name__ == "__main__":
    sys.exit(main())
