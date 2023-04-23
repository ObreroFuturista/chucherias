import requests
import json
import os
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

def download_image(image_url, save_path):
    response = requests.get(image_url)

    with open(save_path, "wb") as f:
        f.write(response.content)

def download_images(taxon_id, save_directory, place_ids, per_page=100, max_pages=None, research_only=False, num_threads=10):
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)

    base_url = "https://api.inaturalist.org/v1/observations"
    page = 1

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        for place_id in place_ids:
            while True:
                if max_pages is not None and page > max_pages:
                    break

                params = {
                    "taxon_id": taxon_id,
                    "per_page": per_page,
                    "page": page,
                    "photos": "true",
                    "place_id": place_id,
                }

                if research_only:
                    params["quality_grade"] = "research"

                response = requests.get(base_url, params=params)
                data = json.loads(response.text)

                if not data["results"]:
                    break

                futures = []

                for result in data["results"]:
                    if research_only and result["quality_grade"] != "research":
                        continue

                    genus = result["taxon"]["name"].split()[0]  # Extrae el nombre del género
                    research_suffix = "_research" if result["quality_grade"] == "research" else ""
                    for photo in result["photos"]:
                        thumbnail_url = photo["url"]
                        image_url = thumbnail_url.replace("square", "large")
                        image_id = photo["id"]

                        save_path = os.path.join(save_directory, f"{genus}_{image_id}{research_suffix}.jpg")
                        futures.append(executor.submit(download_image, image_url, save_path))

                for future in as_completed(futures):
                    try:
                        future.result()
                    except Exception as e:
                        print(f"Error al descargar la imagen: {e}")

                page += 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Descarga imágenes de iNaturalist para un taxón específico.")
    parser.add_argument("taxon_id", type=int, help="ID del taxón para descargar imágenes")
    parser.add_argument("save_directory", help="Directorio donde se guardarán las imágenes")
    parser.add_argument("--place_ids", nargs='+', type=int, help="IDs de los lugares para restringir las observaciones")
    parser.add_argument("--per_page", type=int, default=100, help="Número de observaciones por página (predeterminado: 100)")
    parser.add_argument("--max_pages", type=int, default=None, help="Número máximo de páginas para descargar (predeterminado: todas)")
    parser.add_argument("--research_only", action="store_true", help="Incluir solo observaciones con grado de investigación")
    parser.add_argument("--num_threads", type=int, default=10, help="Número de hilos para paralelizar las descargas (predeterminado: 10)")

    args = parser.parse_args()

    download_images(args.taxon_id, args.save_directory, args.place_ids, per_page=args.per_page, max_pages=args.max_pages, research_only=args.research_only, num_threads=args.num_threads)
