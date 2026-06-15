# -*- coding: utf-8 -*-
# ## ##########################################################################
#
# Authors: Miguel Nahuatlato
# License: MIT
#
# ========================================================
# Calculadora de la Ley de Paschen  —  v1
# ========================================================
#
# Uso:
#    py paschen_calculator.py
#========================================================
# ## ##########################################################################

import math
from datetime import datetime
import os

# ==============================
# CONSTANTES DE IONIZACIÓN
# ==============================
IONIZATION_DATA = {
    'H2': {'A': 5, 'B': 130},
    'N2': {'A': 12, 'B': 342},
    'CO2': {'A': 20, 'B': 466},
    'He': {'A': 3, 'B': 34},
    'Hg': {'A': 20, 'B': 370},
    'air': {'A': 15, 'B': 365},  # aire
}

# Diccionario para almacenar gases personalizados agregados por el usuario
GASES_PERSONALIZADOS = {}

# ==============================
# FUNCIÓN PARA LIMPIAR PANTALLA
# ==============================
def limpiar_pantalla():
    """Limpia la pantalla de la terminal"""
    os.system('cls' if os.name == 'nt' else 'clear')
    print("✓ Pantalla limpiada\n")

# ==============================
# FUNCIÓN PARA NORMALIZAR NOMBRES DE GASES
# ==============================
def normalizar_nombre_gas(nombre):
    """
    Normaliza el nombre del gas para aceptar minúsculas/mayúsculas
    pero mantiene los números exactos (N2 ≠ N)
    """
    # Convertir a mayúsculas solo las letras, mantener números
    nombre_normalizado = ''
    for i, char in enumerate(nombre):
        if char.isalpha():
            # Primera letra siempre mayúscula
            if i == 0 or nombre[i-1].isdigit():
                nombre_normalizado += char.upper()
            else:
                # Letras siguientes minúsculas (para casos como Co, He)
                nombre_normalizado += char.lower()
        else:
            # Mantener números y otros caracteres tal cual
            nombre_normalizado += char
    
    # Casos especiales para asegurar formato correcto
    if nombre_normalizado.upper() == 'HE':
        return 'He'
    elif nombre_normalizado.upper() == 'H2':
        return 'H2'
    elif nombre_normalizado.upper() == 'N2':
        return 'N2'
    elif nombre_normalizado.upper() == 'CO2':
        return 'CO2'
    elif nombre_normalizado.upper() == 'HG':
        return 'Hg'
    elif nombre_normalizado.upper() == 'AIR' or nombre_normalizado.upper() == 'AIRE':
        return 'air'
    
    return nombre_normalizado

def obtener_todos_los_gases():
    """Retorna todos los gases disponibles (predefinidos + personalizados)"""
    todos_gases = {**IONIZATION_DATA, **GASES_PERSONALIZADOS}
    return todos_gases

def agregar_gas_personalizado():
    """Permite al usuario agregar un gas personalizado con sus constantes A y B"""
    print("\n" + "="*60)
    print("AGREGAR GAS PERSONALIZADO")
    print("="*60)
    
    try:
        nombre = input("Nombre del gas (ej: Ar, O2, Ne): ").strip()
        nombre_normalizado = normalizar_nombre_gas(nombre)
        
        # Verificar si ya existe
        if nombre_normalizado in IONIZATION_DATA:
            print(f"\n⚠️  El gas '{nombre_normalizado}' ya existe en la base de datos.")
            return False
        
        if nombre_normalizado in GASES_PERSONALIZADOS:
            print(f"\n⚠️  El gas '{nombre_normalizado}' ya fue agregado anteriormente.")
            sobrescribir = input("¿Deseas sobrescribir sus valores? (s/n): ").strip().lower()
            if sobrescribir != 's':
                return False
        
        A = float(input(f"Constante A para {nombre_normalizado} (1/(Torr·cm)): "))
        B = float(input(f"Constante B para {nombre_normalizado} (V/(Torr·cm)): "))
        
        GASES_PERSONALIZADOS[nombre_normalizado] = {'A': A, 'B': B}
        
        print(f"\n✓ Gas '{nombre_normalizado}' agregado exitosamente!")
        print(f"  A = {A} 1/(Torr·cm)")
        print(f"  B = {B} V/(Torr·cm)")
        
        return True
        
    except ValueError as e:
        print(f"\n❌ Error: {e}")
        return False
    except Exception as e:
        print(f"\n❌ Error inesperado: {e}")
        return False

# ==============================
# FUNCIÓN PARA NORMALIZAR UNIDADES DE PRESIÓN
# ==============================
def normalizar_unidad_presion(unidad):
    """
    Normaliza la unidad de presión para aceptar minúsculas/mayúsculas
    """
    unidad_upper = unidad.upper()
    
    # Mapeo de posibles variaciones a la unidad correcta
    mapeo_unidades = {
        'PA': 'Pa',
        'TORR': 'Torr',
        'ATM': 'atm',
        'MBAR': 'mbar',
        'PSI': 'psi',
        'MMHG': 'mmHg',
        'KG/CM²': 'kg/cm²',
        'KG/CM2': 'kg/cm²',
        'KGCM2': 'kg/cm²',
        'AT': 'kg/cm²',
        'INHG': 'inHg',
        'IN HG': 'inHg',
        'IN_HG': 'inHg'
    }
    
    if unidad_upper in mapeo_unidades:
        return mapeo_unidades[unidad_upper]
    else:
        raise ValueError(f"Unidad '{unidad}' no reconocida. Unidades válidas: Pa, Torr, atm, mbar, psi, mmHg, inHg, kg/cm²")

# ==============================
# FUNCIONES DE CONVERSIÓN DE PRESIÓN
# ==============================
def convertir_presion(valor, unidad_entrada):
    """Convierte presión de cualquier unidad a todas las demás"""
    # Normalizar la unidad de entrada
    unidad_entrada = normalizar_unidad_presion(unidad_entrada)
    
    # Primero convertir todo a Pascal (Pa)
    conversiones_a_pa = {
        'Pa': 1,
        'Torr': 133.322,
        'atm': 101325,
        'mbar': 100,
        'psi': 6894.76,
        'mmHg': 133.322,
        'kg/cm²': 98066.5,
        'inHg': 3386.39
    }
    
    # Convertir a Pascal
    valor_pa = valor * conversiones_a_pa[unidad_entrada]
    
    # Convertir de Pascal a todas las unidades
    resultados = {}
    for unidad, factor in conversiones_a_pa.items():
        resultados[unidad] = valor_pa / factor
    
    return resultados

def mostrar_conversiones(valor, unidad_entrada):
    """Muestra todas las conversiones de presión"""
    resultados = convertir_presion(valor, unidad_entrada)
    
    print(f"\n{'='*60}")
    print(f"Conversiones de Presión: {valor} {unidad_entrada}")
    print(f"{'='*60}")
    print(f"  Pascal (Pa):       {resultados['Pa']:.6f} Pa")
    print(f"  Torr:              {resultados['Torr']:.6f} Torr")
    print(f"  Atmósferas (atm):  {resultados['atm']:.6f} atm")
    print(f"  Milibar (mbar):    {resultados['mbar']:.6f} mbar")
    print(f"  PSI:               {resultados['psi']:.6f} psi")
    print(f"  mmHg:              {resultados['mmHg']:.6f} mmHg")
    print(f"  inHg:              {resultados['inHg']:.6f} inHg")
    print(f"  kg/cm² (at):       {resultados['kg/cm²']:.6f} kg/cm²")
    print(f"{'='*60}\n")
    
    return resultados

# ==============================
# FUNCIONES PARA CÁLCULO DE A Y B
# ==============================
def calcular_A_B_mezcla(porcentajes):
    """
    Calcula los valores A y B de la mezcla usando porcentajes de concentración
    porcentajes: dict {'gas': porcentaje} donde porcentaje está en %
    """
    # Obtener todos los gases disponibles
    todos_gases = obtener_todos_los_gases()
    
    # Convertir porcentajes a fracciones molares
    suma_porcentajes = sum(porcentajes.values())
    
    if not math.isclose(suma_porcentajes, 100.0, rel_tol=1e-3):
        print(f"ADVERTENCIA: Los porcentajes suman {suma_porcentajes:.2f}%, no 100%")
        print("Se normalizarán los porcentajes...")
    
    # Convertir a fracciones molares (dividir entre 100 y normalizar)
    fracciones_molares = {gas: (porc/suma_porcentajes) for gas, porc in porcentajes.items()}
    
    A_mezcla = 0
    B_mezcla = 0
    
    print(f"\n{'='*60}")
    print("Cálculo de constantes A y B de la mezcla")
    print(f"{'='*60}")
    
    for gas, fraccion in fracciones_molares.items():
        if gas not in todos_gases:
            raise ValueError(f"Gas '{gas}' no está en la base de datos")
        
        A_gas = todos_gases[gas]['A']
        B_gas = todos_gases[gas]['B']
        
        A_mezcla += fraccion * A_gas
        B_mezcla += fraccion * B_gas
        
        porcentaje_original = porcentajes[gas]
        porcentaje_normalizado = fraccion * 100
        
        print(f"{gas}:")
        print(f"  Porcentaje ingresado: {porcentaje_original:.2f}%")
        if not math.isclose(suma_porcentajes, 100.0, rel_tol=1e-3):
            print(f"  Porcentaje normalizado: {porcentaje_normalizado:.2f}%")
        print(f"  Fracción molar: {fraccion:.6f}")
        print(f"  A = {A_gas} 1/(Torr·cm), B = {B_gas} V/(Torr·cm)")
        print(f"  Contribución a A: {fraccion * A_gas:.6f}")
        print(f"  Contribución a B: {fraccion * B_gas:.6f}")
        print()
    
    print(f"{'='*60}")
    print(f"RESULTADOS DE LA MEZCLA:")
    print(f"  A (mezcla) = {A_mezcla:.6f} 1/(Torr·cm)")
    print(f"  B (mezcla) = {B_mezcla:.6f} V/(Torr·cm)")
    print(f"{'='*60}\n")
    
    return A_mezcla, B_mezcla, fracciones_molares, porcentajes

# ==============================
# FUNCIONES DE PASCHEN
# ==============================
def paschen_voltage(p, A, B, d, gamma):
    """Calcula el voltaje de ruptura usando la ley de Paschen"""
    numerator = B * p * d
    denominator = math.log(A * p * d) - math.log(math.log(1 + 1/gamma))
    return numerator / denominator

def solve_for_p(Vb_target, A, B, d, gamma, p_min=1e-6, p_max=1e6, tol=1e-6, max_iter=1000):
    """Resuelve para encontrar la presión dado un voltaje de ruptura"""
    def f(p):
        return paschen_voltage(p, A, B, d, gamma) - Vb_target
    
    if f(p_min) * f(p_max) >= 0:
        raise ValueError("No se puede encontrar la solución en el rango dado.")
    
    for _ in range(max_iter):
        p_mid = (p_min + p_max) / 2
        if abs(f(p_mid)) < tol:
            return p_mid
        if f(p_min) * f(p_mid) < 0:
            p_max = p_mid
        else:
            p_min = p_mid
    
    raise RuntimeError("No convergió la solución.")

# ==============================
# FUNCIÓN PARA GENERAR REPORTE
# ==============================
def generar_reporte(tipo_calculo, datos):
    """Genera un archivo de texto con el reporte del cálculo"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    nombre_archivo = f"reporte_paschen_{timestamp}.txt"
    
    # Obtener todos los gases disponibles
    todos_gases = obtener_todos_los_gases()
    
    with open(nombre_archivo, 'w', encoding='utf-8') as f:
        f.write("="*70 + "\n")
        f.write("REPORTE DE CÁLCULOS - LEY DE PASCHEN\n")
        f.write("="*70 + "\n")
        f.write(f"Fecha y hora: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Tipo de cálculo: {tipo_calculo}\n")
        f.write("="*70 + "\n\n")
        
        if tipo_calculo == "Cálculo de A y B por fracciones molares":
            f.write("COMPOSICIÓN DE LA MEZCLA:\n")
            f.write("-"*70 + "\n")
            
            suma_porcentajes = sum(datos['porcentajes'].values())
            
            for gas, porcentaje in datos['porcentajes'].items():
                fraccion = datos['fracciones_molares'][gas]
                f.write(f"{gas}:\n")
                f.write(f"  Porcentaje: {porcentaje:.2f}%\n")
                if not math.isclose(suma_porcentajes, 100.0, rel_tol=1e-3):
                    f.write(f"  Porcentaje normalizado: {fraccion * 100:.2f}%\n")
                f.write(f"  Fracción molar: {fraccion:.6f}\n")
                f.write(f"  A = {todos_gases[gas]['A']} 1/(Torr·cm)\n")
                f.write(f"  B = {todos_gases[gas]['B']} V/(Torr·cm)\n")
                f.write(f"  Contribución a A: {fraccion * todos_gases[gas]['A']:.6f}\n")
                f.write(f"  Contribución a B: {fraccion * todos_gases[gas]['B']:.6f}\n")
                f.write("\n")
            
            if not math.isclose(suma_porcentajes, 100.0, rel_tol=1e-3):
                f.write(f"NOTA: Suma de porcentajes = {suma_porcentajes:.2f}%\n")
                f.write("Los porcentajes fueron normalizados a 100%\n\n")
            
            f.write("="*70 + "\n")
            f.write("RESULTADOS DE LA MEZCLA:\n")
            f.write("-"*70 + "\n")
            f.write(f"A (mezcla) = {datos['A_mezcla']:.6f} 1/(Torr·cm)\n")
            f.write(f"B (mezcla) = {datos['B_mezcla']:.6f} V/(Torr·cm)\n")
            f.write("="*70 + "\n\n")
        
        if tipo_calculo == "Cálculo de presión óptima (Ley de Paschen)":
            f.write("PARÁMETROS DE ENTRADA:\n")
            f.write("-"*70 + "\n")
            f.write(f"Voltaje de ruptura (Vb): {datos['Vb']} V\n")
            f.write(f"Constante A: {datos['A']} 1/(Torr·cm)\n")
            f.write(f"Constante B: {datos['B']} V/(Torr·cm)\n")
            f.write(f"Distancia entre electrodos (d): {datos['d']} cm\n")
            f.write(f"Coeficiente gamma: {datos['gamma']}\n\n")
            
            f.write("="*70 + "\n")
            f.write("RESULTADO:\n")
            f.write("-"*70 + "\n")
            f.write(f"Presión calculada: {datos['presion_torr']:.6f} Torr\n\n")
            
            f.write("CONVERSIONES DE PRESIÓN:\n")
            f.write("-"*70 + "\n")
            conversiones = convertir_presion(datos['presion_torr'], 'Torr')
            f.write(f"  Pascal (Pa):       {conversiones['Pa']:.6f} Pa\n")
            f.write(f"  Torr:              {conversiones['Torr']:.6f} Torr\n")
            f.write(f"  Atmósferas (atm):  {conversiones['atm']:.6f} atm\n")
            f.write(f"  Milibar (mbar):    {conversiones['mbar']:.6f} mbar\n")
            f.write(f"  PSI:               {conversiones['psi']:.6f} psi\n")
            f.write(f"  mmHg:              {conversiones['mmHg']:.6f} mmHg\n")
            f.write(f"  inHg:              {conversiones['inHg']:.6f} inHg\n")
            f.write(f"  kg/cm² (at):       {conversiones['kg/cm²']:.6f} kg/cm²\n")
            f.write("="*70 + "\n\n")
        
        f.write("\n" + "="*70 + "\n")
        f.write("FIN DEL REPORTE\n")
        f.write("="*70 + "\n")
    
    print(f"\n✓ Reporte generado: {nombre_archivo}\n")
    return nombre_archivo

# ==============================
# MENÚ PRINCIPAL
# ==============================
def menu_principal():
    """Muestra el menú principal y maneja la navegación"""
    while True:
        print("\n" + "="*60)
        print("CALCULADORA DE LEY DE PASCHEN")
        print("="*60)
        print("1. Limpiar pantalla")
        print("2. Conversión de unidades de presión")
        print("3. Calcular A y B por porcentajes de concentración")
        print("4. Calcular presión óptima (Ley de Paschen)")
        print("5. Salir")
        print("="*60)
        
        opcion = input("Selecciona una opción (1-5): ").strip()
        
        if opcion == '1':
            limpiar_pantalla()
        elif opcion == '2':
            opcion_conversion_presion()
        elif opcion == '3':
            opcion_calcular_A_B()
        elif opcion == '4':
            opcion_calcular_presion()
        elif opcion == '5':
            print("\n¡Hasta luego!\n")
            break
        else:
            print("\n❌ Opción no válida. Intenta de nuevo.")

def opcion_conversion_presion():
    """Opción 2: Conversión de unidades de presión"""
    while True:
        print("\n" + "="*60)
        print("CONVERSIÓN DE UNIDADES DE PRESIÓN")
        print("="*60)
        print("Unidades disponibles: Pa, Torr, atm, mbar, psi, mmHg, inHg, kg/cm²")
        print("(Puedes escribirlas en mayúsculas o minúsculas)")
        
        try:
            valor = float(input("\nIngresa el valor de presión: "))
            unidad_input = input("Ingresa la unidad: ").strip()
            
            # Normalizar y mostrar conversiones
            mostrar_conversiones(valor, unidad_input)
            
            # Preguntar si desea repetir la conversión
            repetir = input("¿Deseas hacer otra conversión? (s/n): ").strip().lower()
            
            if repetir != 's':
                limpiar_pantalla()
                print("\nVolviendo al menú principal...")
                break
            
            limpiar_pantalla()
            
        except ValueError as e:
            print(f"\n{'='*60}")
            print("❌ ERROR DETECTADO")
            print(f"{'='*60}")
            print(f"Descripción: {e}")
            print(f"{'='*60}")
            
            reintentar = input("\n¿Deseas intentarlo de nuevo? (s/n): ").strip().lower()
            
            if reintentar != 's':
                limpiar_pantalla()
                print("\nVolviendo al menú principal...")
                break
            
            limpiar_pantalla()
            
        except Exception as e:
            print(f"\n{'='*60}")
            print("❌ ERROR INESPERADO")
            print(f"{'='*60}")
            print(f"Descripción: {e}")
            print(f"{'='*60}")
            
            reintentar = input("\n¿Deseas intentarlo de nuevo? (s/n): ").strip().lower()
            
            if reintentar != 's':
                limpiar_pantalla()
                print("\nVolviendo al menú principal...")
                break
            
            limpiar_pantalla()

def opcion_calcular_A_B():
    """Opción 3: Calcular A y B por porcentajes de concentración"""
    while True:
        # Obtener gases disponibles
        todos_gases = obtener_todos_los_gases()
        lista_gases = ', '.join(sorted(todos_gases.keys()))
        
        print("\n" + "="*60)
        print("CÁLCULO DE A Y B POR PORCENTAJES DE CONCENTRACIÓN")
        print("="*60)
        print(f"Gases disponibles: {lista_gases}")
        print("Nota: Ingresa los porcentajes de cada gas (ej: 30, 50, 20)")
        print("      Los porcentajes deberían sumar 100% (se normalizan si no)")
        print("      Puedes escribir los gases en minúsculas (he = He)")
        print("="*60)
        
        try:
            porcentajes = {}
            
            # Pedir número de gases
            n_gases = int(input("\n¿Cuántos gases tiene la mezcla? (1-10): "))
            
            if n_gases < 1 or n_gases > 10:
                raise ValueError("El número de gases debe estar entre 1 y 10")
            
            for i in range(n_gases):
                print(f"\n--- Gas {i+1} ---")
                
                # Primero pedir el porcentaje
                porcentaje = float(input("  Porcentaje (%): "))
                
                if porcentaje < 0:
                    raise ValueError("El porcentaje no puede ser negativo")
                
                # Luego pedir el nombre del gas
                gas_input = input(f"  Nombre del gas ({lista_gases}): ").strip()
                gas = normalizar_nombre_gas(gas_input)
                
                # Verificar si el gas existe
                if gas not in todos_gases:
                    print(f"\n⚠️  El gas '{gas_input}' no está en la base de datos.")
                    agregar = input("¿Deseas agregarlo como gas personalizado? (s/n): ").strip().lower()
                    
                    if agregar == 's':
                        print(f"\nAgregando '{gas}' a la base de datos...")
                        try:
                            A = float(input(f"  Constante A para {gas} (1/(Torr·cm)): "))
                            B = float(input(f"  Constante B para {gas} (V/(Torr·cm)): "))
                            GASES_PERSONALIZADOS[gas] = {'A': A, 'B': B}
                            todos_gases = obtener_todos_los_gases()
                            print(f"  ✓ Gas '{gas}' agregado exitosamente!")
                        except ValueError as e:
                            raise ValueError(f"Error al agregar gas personalizado: {e}")
                    else:
                        raise ValueError(f"Gas '{gas_input}' no reconocido y no fue agregado")
                
                if gas in porcentajes:
                    print(f"  ⚠️  Advertencia: {gas} ya fue ingresado. Se sumará el porcentaje.")
                    porcentajes[gas] += porcentaje
                else:
                    porcentajes[gas] = porcentaje
                
                print(f"  ✓ {gas}: {porcentaje}%")
            
            # Calcular A y B
            A_mezcla, B_mezcla, fracciones_molares, porcentajes_orig = calcular_A_B_mezcla(porcentajes)
            
            # Preguntar si desea generar reporte PRIMERO
            generar = input("\n¿Deseas generar un reporte en .txt? (s/n): ").strip().lower()
            
            if generar == 's':
                datos_reporte = {
                    'porcentajes': porcentajes_orig,
                    'fracciones_molares': fracciones_molares,
                    'A_mezcla': A_mezcla,
                    'B_mezcla': B_mezcla
                }
                generar_reporte("Cálculo de A y B por fracciones molares", datos_reporte)
            
            # LUEGO preguntar si desea repetir
            repetir = input("¿Deseas calcular otra mezcla? (s/n): ").strip().lower()
            
            if repetir != 's':
                limpiar_pantalla()
                print("\nVolviendo al menú principal...")
                break
            
            limpiar_pantalla()
            
        except ValueError as e:
            print(f"\n{'='*60}")
            print("❌ ERROR DETECTADO")
            print(f"{'='*60}")
            print(f"Descripción: {e}")
            print(f"{'='*60}")
            
            reintentar = input("\n¿Deseas intentarlo de nuevo? (s/n): ").strip().lower()
            
            if reintentar != 's':
                limpiar_pantalla()
                print("\nVolviendo al menú principal...")
                break
            
            limpiar_pantalla()
            
        except Exception as e:
            print(f"\n{'='*60}")
            print("❌ ERROR INESPERADO")
            print(f"{'='*60}")
            print(f"Descripción: {e}")
            print(f"{'='*60}")
            
            reintentar = input("\n¿Deseas intentarlo de nuevo? (s/n): ").strip().lower()
            
            if reintentar != 's':
                limpiar_pantalla()
                print("\nVolviendo al menú principal...")
                break
            
            limpiar_pantalla()

def opcion_calcular_presion():
    """Opción 4: Calcular presión óptima usando Ley de Paschen"""
    while True:
        print("\n" + "="*60)
        print("CÁLCULO DE PRESIÓN ÓPTIMA (LEY DE PASCHEN)")
        print("="*60)
        
        try:
            Vb = float(input("\nIngresa el Voltaje de Ruptura Vb (V): "))
            A = float(input("Ingresa la constante A (1/(Torr·cm)): "))
            B = float(input("Ingresa la constante B (V/(Torr·cm)): "))
            d = float(input("Ingresa la distancia entre electrodos d (cm): "))
            gamma = float(input("Ingresa el coeficiente gamma: "))
            
            # Calcular presión
            p = solve_for_p(Vb, A, B, d, gamma)
            
            print(f"\n{'='*60}")
            print(f"RESULTADO:")
            print(f"{'='*60}")
            print(f"Presión calculada: {p:.6f} Torr")
            print(f"{'='*60}\n")
            
            # Mostrar conversiones
            mostrar_conversiones(p, 'Torr')
            
            # Preguntar si desea generar reporte PRIMERO
            generar = input("¿Deseas generar un reporte en .txt? (s/n): ").strip().lower()
            
            if generar == 's':
                datos_reporte = {
                    'Vb': Vb,
                    'A': A,
                    'B': B,
                    'd': d,
                    'gamma': gamma,
                    'presion_torr': p
                }
                generar_reporte("Cálculo de presión óptima (Ley de Paschen)", datos_reporte)
            
            # LUEGO preguntar si desea repetir
            repetir = input("¿Deseas calcular otra presión? (s/n): ").strip().lower()
            
            if repetir != 's':
                limpiar_pantalla()
                print("\nVolviendo al menú principal...")
                break
            
            limpiar_pantalla()
            
        except ValueError as e:
            print(f"\n{'='*60}")
            print("❌ ERROR DETECTADO")
            print(f"{'='*60}")
            print(f"Descripción: {e}")
            print(f"{'='*60}")
            
            reintentar = input("\n¿Deseas intentarlo de nuevo? (s/n): ").strip().lower()
            
            if reintentar != 's':
                limpiar_pantalla()
                print("\nVolviendo al menú principal...")
                break
            
            limpiar_pantalla()
            
        except Exception as e:
            print(f"\n{'='*60}")
            print("❌ ERROR INESPERADO")
            print(f"{'='*60}")
            print(f"Descripción: {e}")
            print(f"{'='*60}")
            
            reintentar = input("\n¿Deseas intentarlo de nuevo? (s/n): ").strip().lower()
            
            if reintentar != 's':
                limpiar_pantalla()
                print("\nVolviendo al menú principal...")
                break
            
            limpiar_pantalla()

# ==============================
# EJECUTAR PROGRAMA
# ==============================
if __name__ == "__main__":
    # Limpiar pantalla al iniciar el programa
    limpiar_pantalla()
    print("Bienvenido a la Calculadora de Ley de Paschen\n")
    menu_principal()