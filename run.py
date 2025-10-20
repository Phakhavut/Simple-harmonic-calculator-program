import os
import subprocess
import sys

def main():
    # ติดตั้ง Flask ถ้ายังไม่มี
    try:
        import flask
    except ImportError:
        print("📦 Installing Flask...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", "requirements.txt"])

    print("\n🚀 Starting Simple Harmonic Calculator...")
    os.system(f"{sys.executable} app.py")

if __name__ == "__main__":
    main()
