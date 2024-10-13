from PIL import Image
import itertools

values = [round(i * 0.05, 2) for i in range(1, 9)] 
combinations = itertools.product(values, repeat=2)



# Lista dei nomi dei file PNG
png_files = [f"C:\\Users\\matte\\Documents\\STO MOD AND SIM\\Progetto Esame\\frame_r1.0_l0.5_k1.0_m0.9_frame_{v1}_{v2}.png" for v1, v2 in combinations]  #frame_r1.0_l0.5_k1.0_m0.4_  #frame_r1.0_l0.5_k1.0_m0.9_ #

# Carica le immagini
images = [Image.open(file) for file in png_files]

# Crea la GIF animata e salvala
images[0].save(
    'C:\\Users\\matte\\Documents\\STO MOD AND SIM\\Progetto Esame\\pdf_r1_l05_k1_m09.gif', 
    save_all=True,  # Salva tutte le immagini
    append_images=images[1:],  # Aggiungi le immagini rimanenti
    duration=500,  # Durata di ogni frame in millisecondi
    loop=0  # 0 significa loop infinito
)



