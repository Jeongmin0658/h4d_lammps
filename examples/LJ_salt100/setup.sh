# --- Define where the "ver0.0" folder is ---
#
search_word="WHERE_YOUR_SOURCES_ARE"
your_ver_folder="YOUR_VER_0.0"
sed -i "s%${search_word}%${your_ver_folder}%" run_gcmc_sample.py
sed -i "s%${search_word}%${your_ver_folder}%" dictionaries.py
