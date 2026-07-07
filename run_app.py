import os
import sys
import streamlit.web.bootstrap as bootstrap

if __name__ == '__main__':
    # Determine base path for frozen executable vs running source script
    if getattr(sys, 'frozen', False):
        base_path = sys._MEIPASS
    else:
        base_path = os.path.dirname(os.path.abspath(__file__))
        
    app_path = os.path.join(base_path, 'App.py')
    
    if not os.path.exists(app_path):
        print(f"Error: Could not find main application script App.py at '{app_path}'")
        sys.exit(1)
        
    # Programmatic Streamlit config flags
    flag_options = {
        "server.port": 8501,
        "global.developmentMode": False,
        "server.headless": False,  # Automatically launches the web browser
        "browser.gatherUsageStats": False,
    }
    
    bootstrap.load_config_options(flag_options=flag_options)
    flag_options["_is_running_with_streamlit"] = True
    bootstrap.run(app_path, "streamlit run", [], flag_options)
